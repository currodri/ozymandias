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
from yt import YTArray, YTQuantity
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
    parser = argparse.ArgumentParser(description='Tracking of galaxy masses across cosmic time.')
    parser.add_argument('model', type=str, nargs='+', help='Simulation names from which extract OZY results.')
    parser.add_argument('--ind', type=int, default=0, help='Index of the galaxy in the ending snapshot.')
    parser.add_argument('--var', type=str, nargs='+', default=['stellar'], help='Galaxy mass types to be plotted.')
    parser.add_argument('--start', type=int, default=1, help='Starting index to look for galaxy.')
    parser.add_argument('--end', type=int, default=1000, help='Ending index to look for galaxy.')
    parser.add_argument('--maxz',type=float,default=2.0, help="Maximum redshift displayed in plot.")
    parser.add_argument('--NUT', type=bool, default=True, help='If True, it looks for NUT as the most massive galaxy (stars) in the last snapshot.')
    args = parser.parse_args()

    if not isinstance(args.model, list):
        args.model = [args.model]

    # Now add all to a simple plot
    fig, axes = plt.subplots(len(args.model), 1, sharex=True, figsize=(5,7), dpi=100, facecolor='w', edgecolor='k')
    for i in range(0, len(args.model)):
        ax = axes[i]
        ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
        ax.set_ylabel(r'$M/M_{\odot}$', fontsize=16)
        ax.set_ylim([8e+6,1e+11])
        ax.set_yscale('log')
        ax.tick_params(labelsize=12)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
        if args.model[i][0] != '/':
            simfolder = os.path.join(os.getcwd(), args.model[i])
        else:
            simfolder = args.model[i]
            args.model[i] = args.model[i].split('/')[-1]
        
        ax.text(0.5, 0.2, args.model[i],transform=ax.transAxes,fontsize=16)
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

        for mass_var in args.var:
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
            cosmo = FlatLambdaCDM(H0=100*sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
            Ob0=sim.simulation.parameters['omega_b'],Tcmb0=2.73)

            # Age of Universe at this redshift
            thubble = cosmo.age(redshift).value
            ugly = False
            for mass_var in args.var:
                m = sim.galaxies[progind].mass[mass_var]
                m = YTQuantity(m, 'code_mass', registry=sim.unit_registry)
                galaxy_masses[mass_var].append(m.in_units('Msun').d)
            galaxy_time.append(thubble)
            bad_mass = np.zeros(len(args.var))
            for k,mass_var in enumerate(args.var):
                if len(galaxy_time) >= 2:
                    if (galaxy_masses[mass_var][-1]-galaxy_masses[mass_var][-2])/galaxy_masses[mass_var][-2] > 20:
                        bad_mass[k] = True
                    else:
                        bad_mass[k] = False
            
            if np.any(bad_mass == True):
                print('Deleting')
                for m_var in args.var:
                    galaxy_masses[m_var][-2] = galaxy_masses[m_var][-1]

            try:
                progind = sim.galaxies[progind].progen_galaxy_star
            except:
                print('I have lost this galaxy in the snapshot %s'%ozyfile)
                progind = -1

        for mass_var in args.var:
            mass = galaxy_masses[mass_var][::-1]
            time = galaxy_time[::-1]
            ax.plot(time, mass, marker='o', markersize=2, label=mass_var)
        if i==0:
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

        if i == 1:    
            ax.legend(loc='best', fontsize=14,frameon=False)
    fig.subplots_adjust(top=0.91, bottom=0.08,right=0.97,hspace=0.0)
    if args.NUT:
        args.ind = 'NUT'
    fig.savefig(os.getcwd()+'/mass_evolution_'+str(args.ind)+'.png', format='png', dpi=200)

