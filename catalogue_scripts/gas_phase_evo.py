"""
This is a simple code that tracks the stellar and gas mass evolution of a particular galaxy across the set of 
snapshots selected.

By: F. Rodriguez Montero (23/02/2021)
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
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
matplotlib.rcParams.update(params)

if __name__ == '__main__':

    # Parse the command line arguments.
    parser = argparse.ArgumentParser(description='Tracking of gas masses in each phase across cosmic time.')
    parser.add_argument('model', type=str, nargs='+', help='Simulation names from which extract OZY results.')
    parser.add_argument('--ind', type=int, default=0, help='Index of the galaxy in the ending snapshot.')
    parser.add_argument('--start', type=int, default=1, help='Starting index to look for galaxy.')
    parser.add_argument('--end', type=int, default=1000, help='Ending index to look for galaxy.')
    parser.add_argument('--maxz',type=float,default=2.0, help="Maximum redshift displayed in plot.")
    parser.add_argument('--NUT', type=bool, default=True, help='If True, it looks for NUT as the most massive galaxy (stars) in the last snapshot.')
    args = parser.parse_args()

    if not isinstance(args.model, list):
        args.model = [args.model]

    gas_vars = ['gas','gas_cold','gas_warm','gas_hot']
    linestyles = [':','--','-','.-']

    # Now add all to a simple plot
    fig, axes = plt.subplots(1, 1, sharex=True, figsize=(7,5), dpi=100, facecolor='w', edgecolor='k')
    for i in range(0, len(args.model)):
        ax = axes
        ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
        ax.set_ylabel(r'$M_{\rm gas,i}/M_{\rm gas,tot}$', fontsize=16)
        ax.set_ylim([0,100])
        ax.tick_params(labelsize=12)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
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
        galaxy_masses = {}
        galaxy_time = []
        galaxy_position = []
        galaxy_radius = []
        halo_radius = []

        for mass_var in gas_vars:
            galaxy_masses[mass_var] = []

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
            # print(sim.galaxies[progind].mass)
            for mass_var in gas_vars:
                m = sim.galaxies[progind].mass[mass_var]
                galaxy_masses[mass_var].append(m.in_units('Msun').d)
            galaxy_time.append(thubble)

            try:
                progind = sim.galaxies[progind].progen_galaxy_star
            except:
                print('I have lost this galaxy in the snapshot %s'%ozyfile)
                progind = -1

        tot_mass = np.array(galaxy_masses['gas'][::-1])
        # COLD
        mass_cold = np.array(galaxy_masses['gas_cold'][::-1])
        print(mass_cold)
        time = np.array(galaxy_time[::-1])
        ax.plot(time, 100*mass_cold/tot_mass, marker='o', markersize=2, label=args.model[i]+' Cold', color='b',linestyle=linestyles[i])

        # WARM
        mass_warm = np.array(galaxy_masses['gas_warm'][::-1])
        ax.plot(time, 100*(mass_cold + mass_warm)/tot_mass, marker='o', markersize=2, label=args.model[i]+' Warm', color='orange',linestyle=linestyles[i])

        # HOT
        mass_hot = np.array(galaxy_masses['gas_hot'][::-1])
        ax.plot(time, 100*(mass_cold + mass_warm + mass_hot)/tot_mass, marker='o', markersize=2, label=args.model[i]+' Hot', color='r',linestyle=linestyles[i])


        # Add top ticks for redshift
        axR = ax.twiny()
        maxt = cosmo.age(args.maxz).value
        ax.set_xlim(0.0, maxt)
        axR.set_xlim(0.0, maxt)
        topticks1 = np.array([2.0, 3.0, 4.0, 6.0])
        topticks1 = topticks1[topticks1 >= args.maxz]
        topticks2 = cosmo.age(topticks1).value
        axR.set_xticklabels(topticks1)
        axR.set_xticks(topticks2)
        axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
        axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
        axR.set_xlabel(r'$z$', fontsize=16)
        axR.tick_params(labelsize=12)


        ax.legend(loc='best', fontsize=10,frameon=False,ncol=len(args.model))
    fig.subplots_adjust(top=0.91, bottom=0.08,right=0.97,hspace=0.0)
    if args.NUT:
        args.ind = 'NUT'
    fig.savefig(os.getcwd()+'/gas_phase_evolution_'+str(args.ind)+'.png', format='png', dpi=200)

