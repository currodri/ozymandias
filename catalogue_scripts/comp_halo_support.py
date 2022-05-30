"""
This code allows the study of the main sources of gas support across cosmic time
in the halo and the comparison between simulations.

By: F. Rodriguez Montero (30/05/2022)
"""

# Import required libraries
from numpy.core.fromnumeric import var
import ozy
import numpy as np
import os
import sys
import subprocess
import argparse
from astropy.cosmology import FlatLambdaCDM, z_at_value
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/amr')
from amr2 import io_ramses,filtering,amr_integrator
from ozy.utils import init_region,init_filter
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
sns.set(style="white")

if __name__ == '__main__':

    # Parse the command line arguments.
    parser = argparse.ArgumentParser(description='Tracking of halo support across cosmic time.')
    parser.add_argument('model', type=str, nargs='+', help='Simulation names from which extract OZY results.')
    parser.add_argument('--ind', type=int, default=0, help='Index of the galaxy in the ending snapshot.')
    parser.add_argument('--zinit',type=float,default=2.0, help="Redshift at which to start the tracking.")
    parser.add_argument('--zend',type=float,default=6.0, help="Redshift at which to end the tracking.")
    parser.add_argument('--NUT', type=bool, default=True, help='If True, it looks for NUT as the most massive galaxy (stars) in the last snapshot.')
    args = parser.parse_args()

    if not isinstance(args.model, list):
        args.model = [args.model]


    # Now add all to a simple plot
    fig, axes = plt.subplots(1, 1, sharex=True, figsize=(6,5), dpi=100, facecolor='w', edgecolor='k')
    ax = axes
    ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
    ax.set_ylabel(r'$ - \nabla P_X /  \rho_{\rm gas}*\nabla \phi $', fontsize=16)
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

        os.chdir(simfolder)
        result = subprocess.run('IDtoZetas.out -ask '+str(float(args.zinit)), shell=True,stdout=subprocess.PIPE)
        indexinit = int(result.stdout.decode('utf-8').split('is')[1].split('(z')[0])
        result = subprocess.run('IDtoZetas.out -ask '+str(float(args.zend)), shell=True,stdout=subprocess.PIPE)
        indexend = int(result.stdout.decode('utf-8').split('is')[1].split('(z')[0])
        os.chdir('../')

        files = os.listdir(groupspath)

        ozyfiles = []

        for f in files:
            if f.startswith('ozy_') and indexinit <= int(f[4:-5]) <= indexend:
                ozyfiles.append(f)
        
        ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
        ther_support = []
        cr_support = []
        galaxy_time = []
        print(ozyfiles, indexinit, indexend)
        progind = args.ind
        if args.NUT:
            sim = ozy.load(os.path.join(groupspath, ozyfiles[0]))
            virial_mass = [i.virial_quantities['mass'] for i in sim.galaxies]
            progind = np.argmax(virial_mass)

        for j in range(0, len(ozyfiles)):
            ozyfile = ozyfiles[j]
            # if ozyfile == 'ozy_00012.hdf5' and args.model[i]=='cosmoNUThd':
            #     continue
            sim = ozy.load(os.path.join(groupspath, ozyfile))

            if progind == -1:
                break
            else:
                gal = sim.galaxies[progind]

            # Initialise simulation parameters
            redshift = sim.simulation.redshift
            h = sim.simulation.hubble_constant
            cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                                    Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

            # Age of Universe at this redshift
            thubble = cosmo.age(redshift).value

            # Compute galaxy and halo metallicities
            output_path = sim.simulation.fullpath
            rvir =  sim.halos[gal.parent_halo_index].virial_quantities['radius']
            selected_reg = init_region(gal,'sphere', rmin=(0.2,'rvir'),rmax=(rvir.to('kpc').d,'kpc'))
            print('rvir,ozyfile: ',rvir,ozyfile)
            all_filter = init_filter('none','all',gal)
            glob_attrs = amr_integrator.amr_region_attrs()
            if sim.simulation.physics['cr']:
                glob_attrs.nvars = 2
            else:
                glob_attrs.nvars = 1
            
            glob_attrs.nwvars = 1
            glob_attrs.nfilter = 1
            amr_integrator.allocate_amr_regions_attrs(glob_attrs)
            if sim.simulation.physics['cr']:
                glob_attrs.varnames.T.view('S128')[0] = 'grav_therpfrsphere'.ljust(128)
                glob_attrs.varnames.T.view('S128')[1] = 'grav_crpfrsphere'.ljust(128)
            else:
                glob_attrs.varnames.T.view('S128')[0] = 'grav_therpfrsphere'.ljust(128)
            glob_attrs.wvarnames.T.view('S128')[0] = 'volume'.ljust(128)
            glob_attrs.filters[0] = all_filter

            # Begin integration
            amr_integrator.integrate_region(output_path,selected_reg,glob_attrs)
            ther_sp = glob_attrs.data[0,0,0,0]
            ther_support.append(ther_sp)
            if sim.simulatiom.physics['cr']:
                cr_sp = glob_attrs.data[1,0,0,0]
                cr_support.append(cr_sp)
            galaxy_time.append(thubble)

            try:
                progind = gal.progen_galaxy_star
            except:
                print('I have lost this galaxy in the snapshot %s'%ozyfile)
                progind = -1

        ther_support = np.asarray(ther_support[::-1])
        
        time = np.asarray(galaxy_time[::-1])
        if not sim.simulation.physics['cr'] and not sim.simulation.physics['magnetic']:
            ax.plot(time,ther_support, marker='s', markeredgecolor='b', markerfacecolor='none',linestyle=':',color='b',label=args.model[i])
        elif not sim.simulation.physics['cr'] and sim.simulation.physics['magnetic']:
            ax.plot(time,ther_support, marker='s', markeredgecolor='m', markerfacecolor='none',linestyle=':',color='m',label=args.model[i])
        elif sim.simulation.physics['cr']:
            cr_support = np.asarray(cr_support[::-1])
            ax.plot(time,cr_support+ther_support, marker='s', markeredgecolor='g', markerfacecolor='g',label=args.model[i],linestyle='-',color='g')
            ax.plot(time,cr_support, marker='s', markeredgecolor='g', markerfacecolor='g',linestyle='-.',color='g')
            ax.plot(time,ther_support, marker='s', markeredgecolor='g', markerfacecolor='none',linestyle=':',color='g')
        if i==0:
            # Add top ticks for redshift
            axR = ax.twiny()
            maxt = cosmo.age(args.zend).value
            ax.set_xlim(0.0, maxt)
            axR.set_xlim(0.0, maxt)
            topticks1 = np.array([1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0])
            topticks1 = topticks1[topticks1 >= args.zend]
            topticks2 = cosmo.age(topticks1).value
            axR.set_xticklabels(topticks1)
            axR.set_xticks(topticks2)
            axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
            axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
            axR.set_xlabel(r'$z$', fontsize=16)
            axR.tick_params(labelsize=12)

    first_legend = ax.legend(loc='upper right', fontsize=14,frameon=False)
    ax.add_artist(first_legend)
    dummy_lines = []

    dummy_lines.append(ax.plot([],[], color="black", ls = '-',label = r'$ - (\nabla P_{\rm ther} + \nabla P_{\rm cr}) /  \rho_{\rm gas}*\nabla \phi $')[0])
    dummy_lines.append(ax.plot([],[], color="black", ls = '-.',label = r'$ - \nabla P_{\rm cr} /  \rho_{\rm gas}*\nabla \phi $')[0])
    dummy_lines.append(ax.plot([],[], color="black", ls = ':',label = r'$ - \nabla P_{\rm ther} /  \rho_{\rm gas}*\nabla \phi $')[0])
    second_legend = ax.legend(handles=dummy_lines, loc='lower left', frameon=False, fontsize=10)         
    ax.add_artist(second_legend)
    fig.subplots_adjust(top=0.89, bottom=0.1,right=0.97,hspace=0.0)
    if args.NUT:
        args.ind = 'NUT'
    fig.savefig(os.getcwd()+'/comp_halo_support_'+str(args.ind)+'.png', format='png', dpi=300)