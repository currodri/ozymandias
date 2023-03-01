"""
This is a simple code that tracks the evolution of the magnetic field within 
the galactic region, allowing the comparison between different weights,
phases and simulations.

By: Francisco Rodriguez Montero (02/02/2023)
"""

# Import required libraries
import ozy
import numpy as np
import os
import argparse
from astropy.cosmology import FlatLambdaCDM, z_at_value
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
sns.set(style="white")
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
hfont = {'fontname':'Helvetica'}
matplotlib.rc('text', usetex = True)
matplotlib.rc('font', **{'family' : "serif"})
# params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
# matplotlib.rcParams.update(params)

line_dict = {'cosmoNUThd':'royalblue',
                'cosmoNUThd\_all\_cr10':'lightseagreen',
                'cosmoNUThd\_all\_cr20':'lightcoral',
                'cosmoNUTmhd':'m',
                'cosmoNUTcrmhd':'g',
                'cosmoNUTcrmhd\_nost':'olive',
                'cosmoNUTcrmhd\_noheat':'darkgoldenrod',
                'cosmoNUTcrmhd\_3e29':'firebrick',
                'cosmoNUTrticrmhd':'firebrick'}

names = {'cosmoNUThd':'HD',
        'cosmoNUThd\_all\_cr10':'HDcr10',
        'cosmoNUThd\_all\_cr20':'HDcr20',
        'cosmoNUTmhd':'MHD',
        'cosmoNUTcrmhd':'CRMHD',
        'cosmoNUTcrmhd\_nost':'nsCRMHD',
        'cosmoNUTcrmhd\_noheat':'nhCRMHD',
        'cosmoNUTcrmhd\_3e29':r'CRMHD $\kappa = 3\times 10^{29}$ cm$^2$/s',
        'cosmoNUTrticrmhd':'RTCRiMHD'
        }

if __name__ == '__main__':

    # Parse the command line arguments.
    parser = argparse.ArgumentParser(description='Tracking of galaxy B field mass across cosmic time.')
    parser.add_argument('model', type=str, nargs='+', help='Simulation names from which extract OZY results.')
    parser.add_argument('--type', type=str, default='basic', help='Plot type routine to use.')
    parser.add_argument('--ind', type=int, default=0, help='Index of the galaxy in the ending snapshot.')
    parser.add_argument('--start', type=int, default=1, help='Starting index to look for galaxy.')
    parser.add_argument('--end', type=int, default=1000, help='Ending index to look for galaxy.')
    parser.add_argument('--vars', type=str, nargs='+', default=['stellar'], help='Mass variables to explore')
    parser.add_argument('--maxz',type=float,default=1.5, help="Maximum redshift displayed in plot.")
    parser.add_argument('--NUT', type=bool, default=True, help='If True, it looks for NUT as the most massive galaxy (stars) in the last snapshot.')
    args = parser.parse_args()

    if not isinstance(args.model, list):
        args.model = [args.model]
        
    # Now add all to a simple plot
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(6,5), dpi=100, facecolor='w', edgecolor='k')
    ax[0].set_xlabel(r'$t$ [Gyr]', fontsize=16)
    ax[0].set_ylabel(r'$\langle \vert \vec{B}\vert\rangle_{< 0.2 R_{\rm vir}}$ [G]', fontsize=16)
    #ax[0].set_ylim([8e+6,1e+11])
    ax[0].set_yscale('log')
    ax[0].tick_params(labelsize=12)
    ax[0].xaxis.set_ticks_position('both')
    ax[0].yaxis.set_ticks_position('both')
    ax[0].minorticks_on()
    ax[0].tick_params(which='both',axis="both",direction="in")
    ax[1].set_xlabel(r'$t$ [Gyr]', fontsize=16)
    ax[1].set_ylabel(r'$\langle \vert \vec{B}\vert\rangle_{> 0.2 R_{\rm vir}}$ [G]', fontsize=16)
    #ax[1].set_ylim([8e+6,1e+11])
    ax[1].set_yscale('log')
    ax[1].tick_params(labelsize=12)
    ax[1].xaxis.set_ticks_position('both')
    ax[1].yaxis.set_ticks_position('both')
    ax[1].minorticks_on()
    ax[1].tick_params(which='both',axis="both",direction="in")

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

        files = os.listdir(groupspath)

        ozyfiles = []

        for f in files:
            if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                ozyfiles.append(f)
        
        ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
        galaxy_masses = []
        galaxy_time = []
        galaxy_field = []
        halo_field = []

        progind = args.ind
        if args.NUT:
            sim = ozy.load(os.path.join(groupspath, ozyfiles[0]))
            virial_mass = [i.virial_quantities['mass'] for i in sim.galaxies]
            progind = np.argmax(virial_mass)
        for ozyfile in ozyfiles:
            if progind == -1:
                continue
            # Load OZY file
            sim = ozy.load(os.path.join(groupspath, ozyfile))

            # Initialise simulation parameters
            redshift = sim.simulation.redshift
            h = sim.simulation.hubble_constant
            cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                                    Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

            # Age of Universe at this redshift
            thubble = cosmo.age(redshift).value

            m = sim.galaxies[progind].mass['gas']
            galaxy_masses.append(m.in_units('Msun').d)
            galaxy_time.append(thubble)
            # Density weighted median
            galaxy_field.append(np.array([sim.galaxies[progind].magnetism['magnetic_magnitude'][1,1].to('G'),
                                 sim.galaxies[progind].magnetism['magnetic_magnitude_cold'][1,1].to('G'),
                                 sim.galaxies[progind].magnetism['magnetic_magnitude_warm'][1,1].to('G'),
                                 sim.galaxies[progind].magnetism['magnetic_magnitude_hot'][1,1].to('G')]))
            # Volume weighted median
            halo_field.append(sim.galaxies[progind].magnetism['halo_magnetic_magnitude'][3,1].to('G'))
            bad_mass = False
            if len(galaxy_time) >= 2:
                if (galaxy_masses[-1]-galaxy_masses[-2])/galaxy_masses[-2] > 0.8:
                    bad_mass = True
                else:
                    bad_mass = False
            
            if bad_mass == True:
                print('Deleting')
                galaxy_masses[-2] = galaxy_masses[-1]
                galaxy_field[-2] = galaxy_field[-1]
                halo_field[-2] = halo_field[-1]

            try:
                progind = sim.galaxies[progind].progen_galaxy_star
            except:
                print('I have lost this galaxy in the snapshot %s'%ozyfile)
                progind = -1

        field = np.asarray(galaxy_field[::-1])
        hfield = halo_field[::-1]
        time = galaxy_time[::-1]
        ax[0].plot(time, field[:,0], marker='o', markersize=2, color=line_dict[args.model[i]])
        ax[0].plot(time, field[:,3], marker='o', markersize=2, color=line_dict[args.model[i]],linestyle='--')
        ax[0].plot(time, field[:,2], marker='o', markersize=2, color=line_dict[args.model[i]],linestyle=':')
        ax[0].plot(time, field[:,1], marker='o', markersize=2, color=line_dict[args.model[i]],linestyle='-.')
        
        ax[1].plot(time, hfield, marker='o', markersize=2, color=line_dict[args.model[i]],label=names[args.model[i]],)
        
        
        
        # Add top ticks for redshift
        axR = ax[0].twiny()
        maxt = cosmo.age(args.maxz).value
        ax[0].set_xlim(0.0, maxt)
        axR.set_xlim(0.0, maxt)
        topticks1 = np.array([1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
        topticks1 = topticks1[topticks1 >= args.maxz]
        topticks2 = cosmo.age(topticks1).value
        axR.set_xticklabels(topticks1)
        axR.set_xticks(topticks2)
        axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
        axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
        axR.set_xlabel(r'$z$', fontsize=16)
        axR.tick_params(labelsize=12)
    
    ax[1].legend(loc='best', fontsize=14,frameon=False)
        
    dummy_lines = []

    dummy_lines.append(ax[0].plot([],[], color="black", ls = '-',label = 'All')[0])
    dummy_lines.append(ax[0].plot([],[], color="black", ls = '--',label = 'Cold')[0])
    dummy_lines.append(ax[0].plot([],[], color="black", ls = ':',label = 'Warm')[0])
    dummy_lines.append(ax[0].plot([],[], color="black", ls = '-.',label = 'Hot')[0])
    second_legend = ax[0].legend(handles=dummy_lines, loc='lower left', frameon=False, fontsize=10)         
    ax[0].add_artist(second_legend)
    
    fig.subplots_adjust(top=0.91, bottom=0.1,right=0.97,hspace=0.0)
    if args.NUT:
        args.ind = 'NUT'
    fig.savefig(os.getcwd()+'/B_magnitude_evolution_'+str(args.ind)+'.png', format='png', dpi=200)