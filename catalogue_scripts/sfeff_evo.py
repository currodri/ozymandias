"""
This is a simple code that tracks the evolution of the star formation
efficiency given by MTT within 
the galactic region, allowing the comparison between different weights,
phases and simulations.

By: Francisco Rodriguez Montero (07/02/2023)
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
    parser = argparse.ArgumentParser(description='Tracking of galaxy SF efficiency mass across cosmic time.')
    parser.add_argument('model', type=str, nargs='+', help='Simulation names from which extract OZY results.')
    parser.add_argument('--type', type=str, default='basic', help='Plot type routine to use.')
    parser.add_argument('--ind', type=int, default=0, help='Index of the galaxy in the ending snapshot.')
    parser.add_argument('--start', type=int, default=1, help='Starting index to look for galaxy.')
    parser.add_argument('--end', type=int, default=1000, help='Ending index to look for galaxy.')
    parser.add_argument('--maxz',type=float,default=1.5, help="Maximum redshift displayed in plot.")
    parser.add_argument('--NUT', type=bool, default=True, help='If True, it looks for NUT as the most massive galaxy (stars) in the last snapshot.')
    args = parser.parse_args()

    if not isinstance(args.model, list):
        args.model = [args.model]
        
    # Now add all to a simple plot
    fig, ax = plt.subplots(1, 1, sharex=True, figsize=(6,5), dpi=100, facecolor='w', edgecolor='k')
    ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
    ax.set_ylabel(r'$\varepsilon_{\rm ff}$', fontsize=16)
    #ax.set_ylim([8e+6,1e+11])
    #ax.set_yscale('log')
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

        files = os.listdir(groupspath)

        ozyfiles = []

        for f in files:
            if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                ozyfiles.append(f)
        
        ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
        galaxy_masses = []
        galaxy_time = []
        eff = []

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
            # Mass weighted median
            if sim.simulation.physics['cr']:
                eff.append(np.array([sim.galaxies[progind].sf_efficiency['eff_FKmag'][1,1],
                                    sim.galaxies[progind].sf_efficiency['eff_FKmagnocr'][1,1]]))
            elif sim.simulation.physics['magnetic']:
                eff.append(sim.galaxies[progind].sf_efficiency['eff_FKmag'][1,1])
            elif sim.simulation.physics['hydro']:
                eff.append(sim.galaxies[progind].sf_efficiency['eff_FK2'][1,1])
            bad_mass = False
            if len(galaxy_time) >= 2:
                if (galaxy_masses[-1]-galaxy_masses[-2])/galaxy_masses[-2] > 0.8:
                    bad_mass = True
                else:
                    bad_mass = False
            
            if bad_mass == True:
                print('Deleting')
                galaxy_masses[-2] = galaxy_masses[-1]
                eff[-2] = eff[-1]

            try:
                progind = sim.galaxies[progind].progen_galaxy_star
            except:
                print('I have lost this galaxy in the snapshot %s'%ozyfile)
                progind = -1

        eff = np.asarray(eff[::-1])
        time = galaxy_time[::-1]
        if sim.simulation.physics['cr']:
            ax.plot(time, eff[:,0], marker='o', markersize=2, color=line_dict[args.model[i]],label=names[args.model[i]])
            ax.plot(time, eff[:,1], marker='o', markersize=2, color=line_dict[args.model[i]],linestyle='--')
        else:
            ax.plot(time, eff[:], marker='o', markersize=2, color=line_dict[args.model[i]],label=names[args.model[i]])
        
        
        # Add top ticks for redshift
        axR = ax.twiny()
        maxt = cosmo.age(args.maxz).value
        ax.set_xlim(0.0, maxt)
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
    
        ax.legend(loc='best', fontsize=14,frameon=False)
        
    
    fig.subplots_adjust(top=0.91, bottom=0.1,right=0.97,hspace=0.0)
    if args.NUT:
        args.ind = 'NUT'
    fig.savefig(os.getcwd()+'/sf_efficiency_'+str(args.ind)+'.png', format='png', dpi=200)