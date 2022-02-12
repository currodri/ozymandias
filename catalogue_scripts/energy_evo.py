"""
This is a simple code that tracks the energetic evolution of a particular galaxy across the set of 
snapshots selected.

By: F. Rodriguez Montero (15/06/2021)
"""

# Import required libraries
import ozy
from ozy.plot_settings import plotting_dictionary
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
    parser = argparse.ArgumentParser(description='Tracking of galaxy energies across cosmic time.')
    parser.add_argument('model', type=str, nargs='+', help='Simulation names from which extract OZY results.')
    parser.add_argument('--ind', type=int, default=0, help='Index of the galaxy in the ending snapshot.')
    parser.add_argument('--var', type=str, nargs='+', default=['total_thermal'], help='Energies to be plotted.')
    parser.add_argument('--start', type=int, default=1, help='Starting index to look for galaxy.')
    parser.add_argument('--end', type=int, default=1000, help='Ending index to look for galaxy.')
    parser.add_argument('--maxz',type=float,default=1.5, help="Maximum redshift displayed in plot.")
    parser.add_argument('--NUT', type=bool, default=True, help='If True, it looks for NUT as the most massive galaxy (stars) in the last snapshot.')
    args = parser.parse_args()

    if not isinstance(args.model, list):
        args.model = [args.model]

    # Now add all to a simple plot
    fig, axes = plt.subplots(len(args.model), 1, sharex=True, figsize=(5,7), dpi=100, facecolor='w', edgecolor='k')
    for i in range(0, len(args.model)):
        ax = axes[i]
        ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
        ax.set_ylabel(r'$\epsilon$ [erg/g]', fontsize=16)
        ax.set_ylim([8e+9,5e+15])
        ax.set_yscale('log')
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
        galaxy_energies = {}
        galaxy_time = []
        galaxy_masses = []

        for e_var in args.var:
            galaxy_energies[e_var] = []

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

            for e_var in args.var:
                try:
                    e = sim.galaxies[progind].energies[e_var]
                    galaxy_energies[e_var].append(e.in_units('erg*g**-1').d)
                except:
                    pass
            galaxy_time.append(thubble)
            m = sim.galaxies[progind].mass['stellar']
            galaxy_masses.append(m)
            bad_mass = False
            if len(galaxy_time) >= 2:
                if (galaxy_masses[-1]-galaxy_masses[-2])/galaxy_masses[-2] > 0.8:
                    bad_mass = True
                else:
                    bad_mass = False
            
            if bad_mass == True:
                print('Deleting')
                for e_var in args.var:
                    try:
                        galaxy_energies[e_var][-2] = galaxy_energies[e_var][-1]
                    except:
                        pass
            try:
                progind = sim.galaxies[progind].progen_galaxy_star
            except:
                progind = -1

        for e_var in args.var:
            energy = galaxy_energies[e_var][::-1]
            time = galaxy_time[::-1]
            if len(energy) == len(time):
                plotting_def = plotting_dictionary[e_var]
                ax.plot(time, energy, marker='o', markersize=2, label=plotting_def['label'])
        if i==0:
            # Add top ticks for redshift
            axR = ax.twiny()
            maxt = cosmo.age(args.maxz).value
            ax.set_xlim(0.0, maxt)
            axR.set_xlim(0.0, maxt)
            topticks1 = np.array([1.5,2.0, 3.0, 4.0, 6.0])
            topticks1 = topticks1[topticks1 >= args.maxz]
            topticks2 = cosmo.age(topticks1).value
            axR.set_xticklabels(topticks1)
            axR.set_xticks(topticks2)
            axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
            axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
            axR.set_xlabel(r'$z$', fontsize=16)
            axR.tick_params(labelsize=12)

        if i == len(args.model)-1:    
            ax.legend(loc='best', fontsize=14,frameon=False)
    fig.subplots_adjust(top=0.91, bottom=0.08,right=0.97,hspace=0.0)
    if args.NUT:
        args.ind = 'NUT'
    fig.savefig(os.getcwd()+'/specific_energy_evolution_'+str(args.ind)+'.png', format='png', dpi=200)

