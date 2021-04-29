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

        simfolder = os.path.join(os.getcwd(), args.model[i])
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
            stellar_mass = [i.mass['baryon'] for i in sim.galaxies]
            progind = np.argmax(stellar_mass)
            print(args.model[i],progind,sim.simulation.redshift)
        printMass = True
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

            for mass_var in args.var:
                m = sim.galaxies[progind].mass[mass_var]
                m = YTQuantity(m, 'code_mass', registry=sim.unit_registry)
                galaxy_masses[mass_var].append(m.in_units('Msun').d)
            galaxy_time.append(thubble)
            if printMass and redshift >= 4.12:
                print(args.model[i],ozyfile, np.log10(galaxy_masses['stellar'][-1]),progind)
                printMass = False
            galaxy_position.append(sim.galaxies[progind].position)
            r = YTQuantity(sim.galaxies[progind].radius['total'], 'code_length', registry=sim.unit_registry)
            galaxy_radius.append(r.in_units('kpc').d)
            r = YTQuantity(sim.halos[sim.galaxies[progind].parent_halo_index].virial_quantities['radius'],
                            'code_length', registry=sim.unit_registry)
            halo_radius.append(r.in_units('kpc').d)
            try:
                progind = sim.galaxies[progind].progen_galaxy_star
            except:
                progind = -1

        for mass_var in args.var:
            mass = galaxy_masses[mass_var][::-1]
            time = galaxy_time[::-1]
            ax.plot(time, mass, marker='o', markersize=2, label=mass_var)
        if i==0:
            # Add top ticks for redshift
            axR = ax.twiny()
            maxt = cosmo.age(2.0).value
            ax.set_xlim(0.0, maxt)
            axR.set_xlim(0.0, maxt)
            topticks1 = np.array([2.0, 3.0, 4.0, 6.0])
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

    # Now add all to a simple plot
    fig, ax = plt.subplots(1, 1, figsize=(8,5), dpi=100, facecolor='w', edgecolor='k')

    ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
    ax.set_ylabel(r'$r$ [kpc]', fontsize=16)
    ax.tick_params(labelsize=12)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.minorticks_on()
    ax.tick_params(which='both',axis="both",direction="in")

    ax.plot(time, galaxy_radius[::-1], label='gas + stars')
    ax.plot(time, halo_radius[::-1], label=r'20\% DM $r_{\rm vir}$')
    
    # Add top ticks for redshift
    axR = ax.twiny()
    maxt = cosmo.age(2.0).value
    ax.set_xlim(0.0, maxt)
    axR.set_xlim(0.0, maxt)
    topticks1 = np.array([2.0, 3.0, 4.0, 6.0])
    topticks2 = cosmo.age(topticks1).value
    axR.set_xticklabels(topticks1)
    axR.set_xticks(topticks2)
    axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
    axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
    axR.set_xlabel(r'$z$', fontsize=16)
    axR.tick_params(labelsize=12)

    ax.legend(loc='best', fontsize=14)

    fig.savefig(simfolder+'/radius_evolution_'+str(args.ind)+'.png', format='png', dpi=200)

    # Make a simple plot to track the position of the galaxy
    fig, ax = plt.subplots(1, 1, figsize=(8,8), dpi=100, facecolor='w', edgecolor='k')
    ax.set_xlabel(r'$x$ [code length]', fontsize=16)
    ax.set_ylabel(r'$y$ [code length]', fontsize=16)

    galaxy_position = np.array(galaxy_position[::-1])
    time = np.array(galaxy_time[::-1])
    ax.set_xlim([0.99*galaxy_position[:,0].min(), 1.01*galaxy_position[:,0].max()])
    ax.set_ylim([0.99*galaxy_position[:,1].min(), 1.01*galaxy_position[:,1].max()])
    points = np.array([galaxy_position[:,0], galaxy_position[:,1]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    norm = plt.Normalize(time.min(), time.max())
    lc = LineCollection(segments, cmap='viridis', norm=norm)
    # Set the values used for colormapping
    lc.set_array(time)
    lc.set_linewidth(2)
    line = ax.add_collection(lc)
    fig.colorbar(line, ax=ax, label=r'$t$ [Gyr]')

    fig.savefig(simfolder+'/position_track_'+str(args.ind)+'.png', format='png', dpi=200)

