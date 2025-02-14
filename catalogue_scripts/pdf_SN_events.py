"""
This code allows the reading of SN event catalogues from a RAMSES simulation.
If the catalogue is not present, a file needs to be provided in which the list of 
log files to read is given.

By: Curro Rodriguez Montero (currodri@gmail.com)

"""

# Import required libraries
import ozy
from ozy.utils import sn_data_hdf5
from ozy.variables_settings import plotting_dictionary
from unyt import mh,kb
import numpy as np
import os
import argparse
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
sns.set(style="white")
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "Computer Modern Roman",
})
import scipy.stats as st


def contour_hist(x, y, ax, ax_histx, ax_histy, label):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # Peform the kernel density estimate
    xmin,xmax = np.min(x), np.max(x)
    ymin,ymax = np.min(y), np.max(y)
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    rand_indx = np.random.randint(0, x.shape[0],1000)
    values = np.vstack([x[rand_indx], y[rand_indx]])
    print('Using %s data points for the KDE...'%(values.shape[1]))
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    CS = ax.contour(10**xx,10**yy,f, colors=line_dict[label],
                    levels=[1e-2,4e-2,7e-2,1e-1,3e-1],
                    linewidths=[1.5,1.7,1.5,1.7,2.0],
                    linestyles=[':',':','-.','-.','-'])
    ax.plot([],[],color=line_dict[args.model[i]],label=sim2name[label])
    ax.clabel(CS, fontsize=10, inline=True)

    # Histograms
    histx, xedges = np.histogram(x, bins=200)
    ax_histx.stairs(histx/len(x), 10**xedges, color=line_dict[label])
    
    histy, yedges = np.histogram(y, bins=200)
    ax_histy.stairs(histy/len(y), 10**yedges, color=line_dict[label], orientation='horizontal')

if __name__ == '__main__':

    # Parse the command line arguments.
    parser = argparse.ArgumentParser(description='PDF of SN events from RAMSES log files:')
    parser.add_argument('model', type=str, nargs='+', help='Simulation names from which extract SN log.')
    parser.add_argument('--type', type=str, default='1D', help='Plot type: 1D or 2D.')
    parser.add_argument('--flist', type=str, default='log_files.txt', help='Filename of the plain text containing all log files desired.')
    parser.add_argument('--zstart',type=float,default=[10],nargs='+', help="Maximum redshift for the SN event.")
    parser.add_argument('--zend',type=float,default=[1],nargs='+', help="Minimum redshift for the SN event.")
    parser.add_argument('--var', type=str, default='density',help='Which variable to plot the PDF.')
    args = parser.parse_args()

    line_dict = {'cosmoNUThd':'b','cosmoNUTmhd':'m','cosmoNUTcrmhd':'g','cosmoNUTcrmhd\_nost':'olive','cosmoNUTcrmhd\_noheat':'darkkhaki'}
    line_styles = {'cosmoNUThd':':','cosmoNUTmhd':'--','cosmoNUTcrmhd':'-','cosmoNUTcrmhd\_nost':'--','cosmoNUTcrmhd\_noheat':':'}
    sim2name = {'cosmoNUThd':'HD','cosmoNUTmhd':'MHD','cosmoNUTcrmhd':'CRMHD','cosmoNUTcrmhd\_nost':'nostCRMHD','cosmoNUTcrmhd\_noheat':'noheatCRMHD'}
    hatch = {'cosmoNUThd':'-','cosmoNUTmhd':'x','cosmoNUTcrmhd':'*'}

    if not isinstance(args.model, list):
        args.model = [args.model]

    if args.type == '1D':
        # Setup figure
        fig, axes = plt.subplots(1, 1, sharex=True, figsize=(6,3), dpi=100, facecolor='w', edgecolor='k')
        ax = axes
        plot_def = plotting_dictionary[args.var]
        ax.set_xlabel(plot_def['label'], fontsize=16)
        ax.set_ylabel(r'PDF', fontsize=20)
        ax.set_xlim([plot_def['bin_min'],plot_def['bin_max']])
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.tick_params(labelsize=16)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.set_ylim([7e-4,0.1])
        ax.tick_params(which='both',axis="both",direction="in")
        


        for i in range(0, len(args.model)):
            
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
                args.model[i] = args.model[i].replace('_','\_')
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
                args.model[i] = args.model[i].replace('_','\_')
            print(args.model[i])
            
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            
            # Read the file with the names of all log files desired
            logfiles = np.loadtxt(os.path.join(simfolder, args.flist), dtype=str)
            if logfiles.size == 1:
                logfiles = np.array([logfiles])

            # Figure out if we have CRs
            ozy_files = list(filter(lambda x:x[0:3]=='ozy', os.listdir(groupspath)))
            sim = ozy.load(os.path.join(groupspath,ozy_files[0]))

            if sim.simulation.physics['cr']:
                have_crs = True
            else:
                have_crs = False

            # Extract data if it's necessary or update present HDF5 catalogue
            present_dir = os.getcwd()
            os.chdir(simfolder)
            hd = sn_data_hdf5(logfiles,have_crs=have_crs)
            os.chdir(present_dir)

            # Select the desired data given the redshift range given

            z = np.asarray(hd['z'])
            print(z.min(),z.max(), z[:10],z[-10:])
            high_z_ind = int(z[::-1].searchsorted(args.zstart))
            low_z_ind = int(z[::-1].searchsorted(args.zend))
            orig_size = hd['z'].len()
            print('%i SN events between z=%.2f and z=%.2f'%(high_z_ind-low_z_ind, args.zstart, args.zend))

            # And now add this to plot
            if args.var == 'thermal_pressure':
                density = hd['density'][orig_size-high_z_ind:orig_size-low_z_ind]
                temperature = hd['temperature'][orig_size-high_z_ind:orig_size-low_z_ind]
                plotdata = density + temperature + np.log10(kb.to('erg/K'))
            elif args.var == 'density':
                plotdata = hd[args.var][orig_size-high_z_ind:orig_size-low_z_ind] + np.log10(mh.to('g'))
            histy, yedges = np.histogram(plotdata, bins=200)
            ax.stairs(histy/len(plotdata),edges=10**yedges, color=line_dict[args.model[i]]
                      ,label=sim2name[args.model[i]]
                      ,linewidth=1.5)
            hd.close()
            del z

        ax.legend(loc='upper left', fontsize=14,frameon=False,ncol=len(args.model))
        fig.subplots_adjust(top=0.96, bottom=0.22,right=0.98,left=0.15)
        fig.savefig(os.getcwd()+'/pdf_SNevents_'+str(args.var)+'_'+str(args.zstart)+'_'+str(args.zend)+'.pdf', format='pdf', dpi=300)
    
    elif args.type=='2D':

        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        spacing = 0.005


        rect_contour = [left, bottom, width, height]
        rect_histx = [left, bottom + height + spacing, width, 0.2]
        rect_histy = [left + width + spacing, bottom, 0.2, height]

        # And now add this to plot
        xvar = args.var.split('-')[0]
        yvar = args.var.split('-')[1]

        # start with a square Figure
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_axes(rect_contour)
        ax_histx = fig.add_axes(rect_histx, sharex=ax)
        ax_histy = fig.add_axes(rect_histy, sharey=ax)

        plotx = plotting_dictionary[xvar]
        ax.set_xlabel(plotx['label'], fontsize=20)
        ax.set_xscale('log')
        ax.set_xlim([1e-30,8e-20])
        
        ploty = plotting_dictionary[yvar]
        ax.set_ylabel(ploty['label'], fontsize=20)
        ax.set_yscale('log')
        ax.set_ylim([20,8e+8])

        ax.tick_params(labelsize=16)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")

        ax_histx.spines['top'].set_visible(False)
        ax_histx.spines['right'].set_visible(False)
        ax_histx.tick_params(labelsize=16)
        ax_histx.yaxis.set_ticks_position('left')
        ax_histx.xaxis.set_ticks_position('bottom')
        ax_histx.minorticks_on()
        ax_histx.tick_params(which='both',axis="both",direction="in")
        ax_histx.set_yscale('log')
        ax_histx.set_ylim([7e-4,0.1])

        ax_histy.spines['top'].set_visible(False)
        ax_histy.spines['right'].set_visible(False)
        ax_histy.tick_params(labelsize=16)
        ax_histy.yaxis.set_ticks_position('left')
        ax_histy.xaxis.set_ticks_position('bottom')
        ax_histy.minorticks_on()
        ax_histy.tick_params(which='both',axis="both",direction="in")
        ax_histy.set_xscale('log')
        ax_histy.set_xlim([7e-4,0.1])

        # Loop over models
        for i in range(0, len(args.model)):
            
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
                args.model[i] = args.model[i].replace('_','\_')
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
                args.model[i] = args.model[i].replace('_','\_')
            print(args.model[i])
            
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            
            # Read the file with the names of all log files desired
            logfiles = np.loadtxt(os.path.join(simfolder, args.flist), dtype=str)
            if logfiles.size == 1:
                logfiles = np.array([logfiles])

            # Figure out if we have CRs
            ozy_files = list(filter(lambda x:x[0:3]=='ozy', os.listdir(groupspath)))
            sim = ozy.load(os.path.join(groupspath,ozy_files[0]))

            if sim.simulation.physics['cr']:
                have_crs = True
            else:
                have_crs = False

            # Extract data if it's necessary or update present HDF5 catalogue
            present_dir = os.getcwd()
            os.chdir(simfolder)
            hd = sn_data_hdf5(logfiles,have_crs=have_crs)
            os.chdir(present_dir)

            # Select the desired data given the redshift range given

            z = np.asarray(hd['z'])
            print(z.min(),z.max(), z[:10],z[-10:])
            high_z_ind = int(z[::-1].searchsorted(args.zstart))
            low_z_ind = int(z[::-1].searchsorted(args.zend))
            orig_size = hd['z'].len()
            print('%i SN events between z=%.2f and z=%.2f'%(high_z_ind-low_z_ind, args.zstart, args.zend))
            del z

            if xvar == 'thermal_pressure':
                density = hd['density'][orig_size-high_z_ind:orig_size-low_z_ind]
                temperature = hd['temperature'][orig_size-high_z_ind:orig_size-low_z_ind]
                plotdata = density + temperature + np.log10(kb.to('erg/K'))
            elif xvar == 'density':
                plotdata = hd[xvar][orig_size-high_z_ind:orig_size-low_z_ind] + np.log10(mh.to('g'))
            elif xvar == 'temperature':
                plotdata = hd[xvar][orig_size-high_z_ind:orig_size-low_z_ind]
            x = np.copy(plotdata[:])

            if yvar == 'thermal_pressure':
                density = hd['density'][orig_size-high_z_ind:orig_size-low_z_ind]
                temperature = hd['temperature'][orig_size-high_z_ind:orig_size-low_z_ind]
                plotdata = density + temperature + np.log10(kb.to('erg/K'))
            elif yvar == 'density':
                plotdata = hd[yvar][orig_size-high_z_ind:orig_size-low_z_ind] + np.log10(mh.to('g'))
            elif yvar == 'temperature':
                plotdata = hd[yvar][orig_size-high_z_ind:orig_size-low_z_ind]
            y = np.copy(plotdata[:])
            hd.close()
            contour_hist(x, y, ax, ax_histx, ax_histy, args.model[i])
        
        ax.legend(loc='lower left', fontsize=14,frameon=False)
        fig.savefig(os.getcwd()+'/contours_SNevents_'+str(args.var)+'_'+str(args.zstart)+'_'+str(args.zend)+'.png', format='png', dpi=200)
        print('Done')

    elif args.type == '1D_comp':
        # Setup figure
        fig, axes = plt.subplots(2, 1, sharex=True, figsize=(6,6), dpi=300, facecolor='w', edgecolor='k')
        
        ax = axes[0]
        plot_def = plotting_dictionary[args.var]
        ax.set_xlabel(plot_def['label'], fontsize=16)
        ax.set_ylabel(r'PDF', fontsize=20)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.tick_params(labelsize=16)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.set_ylim([5e-4,0.1])
        ax.tick_params(which='both',axis="both",direction="in")
        ax.text(0.05, 0.7, r'$%s<z<%s$'%(args.zend[0],args.zstart[0]),
                            transform=ax.transAxes, fontsize=16,verticalalignment='top',
                            color='black')
        ax = axes[1]
        plot_def = plotting_dictionary[args.var]
        ax.set_xlabel(plot_def['label'], fontsize=16)
        ax.set_ylabel(r'PDF', fontsize=20)
        ax.set_xlim([plot_def['bin_min'],plot_def['bin_max']])
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.tick_params(labelsize=16)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.set_ylim([5e-4,0.1])
        ax.tick_params(which='both',axis="both",direction="in")
        ax.text(0.05, 0.7, r'$%s<z<%s$'%(args.zend[1],args.zstart[1]),
                            transform=ax.transAxes, fontsize=16,verticalalignment='top',
                            color='black')

        for i in range(0, len(args.model)):
            
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
                args.model[i] = args.model[i].replace('_','\_')
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
                args.model[i] = args.model[i].replace('_','\_')
            print(args.model[i])
            
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            
            # Read the file with the names of all log files desired
            logfiles = np.loadtxt(os.path.join(simfolder, args.flist), dtype=str)
            if logfiles.size == 1:
                logfiles = np.array([logfiles])

            # Figure out if we have CRs
            ozy_files = list(filter(lambda x:x[0:3]=='ozy', os.listdir(groupspath)))
            sim = ozy.load(os.path.join(groupspath,ozy_files[0]))

            if sim.simulation.physics['cr']:
                have_crs = True
            else:
                have_crs = False

            # Extract data if it's necessary or update present HDF5 catalogue
            present_dir = os.getcwd()
            os.chdir(simfolder)
            hd = sn_data_hdf5(logfiles,have_crs=have_crs)
            os.chdir(present_dir)

            # Select the desired data given the redshift range given
            z = np.asarray(hd['z'])
            for j in range(0, len(args.zstart)):
            
                high_z_ind = int(z[::-1].searchsorted(args.zstart[j]))
                low_z_ind = int(z[::-1].searchsorted(args.zend[j]))
                orig_size = hd['z'].len()
                print('%i SN events between z=%.2f and z=%.2f'%(high_z_ind-low_z_ind, args.zstart[j], args.zend[j]))

                # And now add this to plot
                if args.var == 'thermal_pressure':
                    density = hd['density'][orig_size-high_z_ind:orig_size-low_z_ind]
                    temperature = hd['temperature'][orig_size-high_z_ind:orig_size-low_z_ind]
                    plotdata = density + temperature + np.log10(kb.to('erg/K'))
                elif args.var == 'density':
                    plotdata = hd[args.var][orig_size-high_z_ind:orig_size-low_z_ind] + np.log10(mh.to('g'))
                histy, yedges = np.histogram(plotdata, bins=200)
                axes[j].stairs(histy/len(plotdata),edges=10**yedges, color=line_dict[args.model[i]]
                        ,label=sim2name[args.model[i]]
                        ,linewidth=1.5)
            hd.close()
            del z

        axes[0].legend(loc='upper left', fontsize=14,frameon=False,ncol=len(args.model))
        fig.subplots_adjust(top=0.96, bottom=0.1,right=0.98,left=0.15,hspace=0,wspace=0)
        fig.savefig(os.getcwd()+'/pdf_SNevents_'+str(args.var)+'_comp.pdf', format='pdf', dpi=300)
    



        



        