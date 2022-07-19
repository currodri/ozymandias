"""
This code allows the comparison of the halo and galaxy metallicity between different
simulations and redshifts.

By: F. Rodriguez Montero (15/02/2022)
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
from ozy.utils import interp_nans
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
sns.set(style="white")
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# hfont = {'fontname':'Helvetica'}
# matplotlib.rc('text', usetex = True)
# matplotlib.rc('font', **{'family' : "serif"})
# params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
# matplotlib.rcParams.update(params)


if __name__ == '__main__':

    # Parse the command line arguments.
    parser = argparse.ArgumentParser(description='Tracking of galaxy metallicity across cosmic time.')
    parser.add_argument('model', type=str, nargs='+', help='Simulation names from which extract OZY results.')
    parser.add_argument('--ind', type=int, default=0, help='Index of the galaxy in the ending snapshot.')
    parser.add_argument('--zinit',type=float,default=2.0, help="Redshift at which to start the tracking.")
    parser.add_argument('--zend',type=float,default=6.0, help="Redshift at which to end the tracking.")
    parser.add_argument('--var', type=str, default='total_thermal', help='Energies to be plotted.')
    parser.add_argument('--NUT', type=bool, default=True, help='If True, it looks for NUT as the most massive galaxy (stars) in the last snapshot.')
    args = parser.parse_args()

    if not isinstance(args.model, list):
        args.model = [args.model]


    # Now add all to a simple plot
    fig, axes = plt.subplots(1, 1, sharex=True, figsize=(7,5), dpi=100, facecolor='w', edgecolor='k')
    ax = axes
    ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
    ax.set_ylabel(r'$u_{\rm CR}$/SFR [erg yr /(cm$^{3}$ $M_{\odot}$]', fontsize=16)
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
        umass = []
        uvolume = []
        ucounts = []
        sfr = []
        galaxy_time = []
        print(ozyfiles, indexinit, indexend)
        progind = args.ind
        if args.NUT:
            sim = ozy.load(os.path.join(groupspath, ozyfiles[0]))
            virial_mass = [i.virial_quantities['mass'] for i in sim.galaxies]
            progind = np.argmax(virial_mass)
            redshift = sim.simulation.redshift
            h = sim.simulation.hubble_constant
            cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                                    Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)
        if not os.path.exists('uvolume.npy'):
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
                sfr.append(gal.sfr['100Myr'].in_units('Msun/yr'))

                # Compute galaxy and halo metallicities
                output_path = sim.simulation.fullpath
                rvir =  sim.halos[gal.parent_halo_index].virial_quantities['radius']
                selected_reg = init_region(gal,'sphere', rmin=(0.0,'rvir'),rmax=(rvir.to('kpc').d,'kpc'))
                print('rvir,ozyfile: ',rvir,ozyfile)
                inner_filter = init_filter('r_sphere/</%.4f/kpc'%(0.2*rvir.to('kpc').d),'galaxy',gal)
                glob_attrs = amr_integrator.amr_region_attrs()
                glob_attrs.nvars = 1
                glob_attrs.nwvars = 3
                glob_attrs.nfilter = 1
                amr_integrator.allocate_amr_regions_attrs(glob_attrs)
                glob_attrs.varnames.T.view('S128')[0] = 'cr_energy_density'.ljust(128)
                glob_attrs.wvarnames.T.view('S128')[0] = 'volume'.ljust(128)
                glob_attrs.wvarnames.T.view('S128')[1] = 'mass'.ljust(128)
                glob_attrs.wvarnames.T.view('S128')[2] = 'counts'.ljust(128)
                glob_attrs.filters[0] = inner_filter

                # Begin integration
                amr_integrator.integrate_region(output_path,selected_reg,glob_attrs)
                uvol = sim.quantity(glob_attrs.data[0,0,0,0],'code_energy_density')
                uvolume.append(uvol.to('erg/cm**3'))
                um = sim.quantity(glob_attrs.data[0,0,1,0],'code_energy_density')
                umass.append(um.to('erg/cm**3'))
                uc = sim.quantity(glob_attrs.data[0,0,2,0],'code_energy_density')
                ucounts.append(uc.to('erg/cm**3'))
                galaxy_time.append(thubble)

                try:
                    progind = gal.progen_galaxy_star
                except:
                    print('I have lost this galaxy in the snapshot %s'%ozyfile)
                    progind = -1

            uvolume = np.asarray(uvolume[::-1])
            umass = np.asarray(umass[::-1])
            ucounts = np.asarray(ucounts[::-1])
            time = np.asarray(galaxy_time[::-1])
            sfr = np.asarray(sfr[::-1])
        if not os.path.exists('uvolume.npy'):
            np.save('uvolume.npy',uvolume)
            np.save('umass.npy',umass)
            np.save('ucounts.npy',ucounts)
            np.save('time.npy',time)
            np.save('sfr.npy',sfr)
        else:
            uvolume = np.load('uvolume.npy')
            umass = np.load('umass.npy')
            ucounts = np.load('ucounts.npy')
            time = np.load('time.npy')
            sfr = np.load('sfr.npy')
        print(np.log10(sfr))
        logsfr = interp_nans(np.log10(sfr))
        sfr = 10**logsfr
        print(uvolume, umass, ucounts, sfr)
        if not sim.simulation.physics['cr'] and not sim.simulation.physics['magnetic']:
            ax.plot(time,uvolume, marker='s', markeredgecolor='b', markerfacecolor='none',linestyle='-',color='b')
        elif not sim.simulation.physics['cr'] and sim.simulation.physics['magnetic']:
            ax.plot(time,uvolume, marker='s', markeredgecolor='m', markerfacecolor='none',linestyle='-',color='m')
        elif sim.simulation.physics['cr']:
            ax.plot(time,uvolume/sfr, marker='s', markeredgecolor='g', markerfacecolor='none',linestyle='-',color='g')
            ax.plot(time,umass/sfr, marker='s', markeredgecolor='g', markerfacecolor='none',linestyle=':',color='g')
            ax.plot(time,ucounts/sfr, marker='s', markeredgecolor='g', markerfacecolor='none',linestyle='-.',color='g')
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

            init_val = ucounts[1]/sfr[1]
            A = init_val * (1+13)**-3
            z = np.linspace(13,5,100)
            model = A*(1+z)**3
            t = cosmo.age(z).value
            print(model,init_val, z, t)
            ax.plot(t, model, 'k:',alpha=0.6)

    dummy_lines = []

    dummy_lines.append(ax.plot([],[], color="black", ls = '-',label = 'Volume-weigthed average')[0])
    dummy_lines.append(ax.plot([],[], color="black", ls = ':',label = 'Mass-weighted average')[0])
    dummy_lines.append(ax.plot([],[], color="black", ls = '-.',label = 'Average')[0])
    second_legend = ax.legend(handles=dummy_lines, loc='lower left', frameon=False, fontsize=10)         
    ax.add_artist(second_legend)
    fig.subplots_adjust(top=0.89, bottom=0.1,right=0.98,left=0.15)
    if args.NUT:
        args.ind = 'NUT'
    fig.savefig(os.getcwd()+'/comp_energy_'+str(args.ind)+'.png', format='png', dpi=300)