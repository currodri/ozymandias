"""
This code allows the reading of SN event catalogues from a RAMSES simulation.
If the catalogue is not present, a file needs to be provided in which the list of 
log files to read is given.

By: Curro Rodriguez Montero (currodri@gmail.com)

"""

# Import required libraries
import ozy
from ozy.utils import sn_data_hdf5
from ozy.plot_settings import plotting_dictionary
from unyt import mh,kb
import numpy as np
import os
import argparse
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
sns.set(style="white")

if __name__ == '__main__':

    # Parse the command line arguments.
    parser = argparse.ArgumentParser(description='PDF of SN events from RAMSES log files:')
    parser.add_argument('model', type=str, nargs='+', help='Simulation names from which extract SN log.')
    parser.add_argument('--type', type=str, default='1D', help='Plot type: 1D or 2D.')
    parser.add_argument('--flist', type=str, default='log_files.txt', help='Filename of the plain text containing all log files desired.')
    parser.add_argument('--zstart',type=float,default=10, help="Maximum redshift for the SN event.")
    parser.add_argument('--zend',type=float,default=1, help="Minimum redshift for the SN event.")
    parser.add_argument('--var', type=str, default='density',help='Which variable to plot the PDF.')
    args = parser.parse_args()

    line_dict = {'cosmoNUThd':'b','cosmoNUTmhd':'m','cosmoNUTcrmhd':'g','cosmoNUTcrmhd\_nost':'olive','cosmoNUTcrmhd\_noheat':'darkkhaki'}
    line_styles = {'cosmoNUThd':':','cosmoNUTmhd':'--','cosmoNUTcrmhd':'-','cosmoNUTcrmhd\_nost':'--','cosmoNUTcrmhd\_noheat':':'}

    if not isinstance(args.model, list):
        args.model = [args.model]

    if args.type == '1D':
        # Setup figure
        fig, axes = plt.subplots(1, 1, sharex=True, figsize=(6,3), dpi=100, facecolor='w', edgecolor='k')
        ax = axes
        plot_def = plotting_dictionary[args.var]
        ax.set_xlabel(plot_def['label_log'], fontsize=16)
        ax.set_ylabel(r'PDF', fontsize=16)
        ax.set_xlim([np.log10(plot_def['vmin']),np.log10(plot_def['vmax'])])
        ax.set_yscale('log')
        ax.tick_params(labelsize=12)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
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
            n,bins, patches = ax.hist(plotdata, bins=100, 
                    range=(np.log10(plot_def['vmin']),np.log10(plot_def['vmax'])),
                    density=True, label=args.model[i], color=line_dict[args.model[i]],
                    histtype='step')
            hd.close()

        ax.legend(loc='lower center', fontsize=14,frameon=False)
        fig.subplots_adjust(top=0.98, bottom=0.22,right=0.98,left=0.13)
        fig.savefig(os.getcwd()+'/pdf_SNevents_'+str(args.var)+'_'+str(args.zstart)+'_'+str(args.zend)+'.png', format='png', dpi=200)
    
    elif args.type=='2D':

        name2label = {'density':r'$\log nH$ [cm$^{-3}$]', 'temperature':r'$\log T$ [K]'}
        # Initialise global lists
        x,y,hue = np.array([]),np.array([]),np.array([])

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

            # And now add this to plot
            xvar = args.var.split('-')[0]
            yvar = args.var.split('-')[1]

            if xvar == 'thermal_pressure':
                density = hd['density'][orig_size-high_z_ind:orig_size-low_z_ind]
                temperature = hd['temperature'][orig_size-high_z_ind:orig_size-low_z_ind]
                plotdata = density + temperature + np.log10(kb.to('erg/K'))
            elif xvar == 'density':
                plotdata = hd[xvar][orig_size-high_z_ind:orig_size-low_z_ind] + np.log10(mh.to('g'))
            elif xvar == 'temperature':
                plotdata = hd[xvar][orig_size-high_z_ind:orig_size-low_z_ind]
            x = np.concatenate((x,plotdata[:]))

            if yvar == 'thermal_pressure':
                density = hd['density'][orig_size-high_z_ind:orig_size-low_z_ind]
                temperature = hd['temperature'][orig_size-high_z_ind:orig_size-low_z_ind]
                plotdata = density + temperature + np.log10(kb.to('erg/K'))
            elif yvar == 'density':
                plotdata = hd[yvar][orig_size-high_z_ind:orig_size-low_z_ind] + np.log10(mh.to('g'))
            elif yvar == 'temperature':
                plotdata = hd[yvar][orig_size-high_z_ind:orig_size-low_z_ind]
            y = np.concatenate((y,plotdata[:]))
            names = len(plotdata[:])*[args.model[i]]
            hue = np.concatenate((hue, np.asarray(names)))

            hd.close()
        # Now construct data frame
        df = pd.DataFrame({name2label[xvar]:x,name2label[yvar]:y,"Model":hue})

        sns.jointplot(data=df, x=name2label[xvar], y=name2label[yvar], hue="Model", kind="kde", marginal_ticks=True, palette=line_dict)
        plt.savefig(os.getcwd()+'/contours_SNevents_'+str(args.var)+'_'+str(args.zstart)+'_'+str(args.zend)+'.png', format='png', dpi=200)


       



        



        