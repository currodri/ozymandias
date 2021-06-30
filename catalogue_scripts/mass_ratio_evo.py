"""
This is a simple code that tracks the DM, stellar and gas mass diference between a set of
runs and a chosen fiducial one (e.g. compare MHD and CRMHD to HD)

By: F. Rodriguez Montero (23/02/2021)
"""
# Import required libraries
import ozy
import numpy as np
import os
import argparse
import itertools
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
    parser = argparse.ArgumentParser(description='Tracking the change in masses due to different models.')
    parser.add_argument('model', type=str, help='Fiducial model used as control.')
    parser.add_argument('--comp', type=str, nargs='+', help='Models used for the comparison to control model.')
    parser.add_argument('--bin',type=float,default=100, help='Size of bins (in Myr) for the moving average (def: 100 Myr).')
    parser.add_argument('--ind', type=int, default=0, help='Index of the galaxy in the ending snapshot.')
    parser.add_argument('--var', type=str, nargs='+', default=['stellar'], help='Galaxy mass types to be plotted.')
    parser.add_argument('--start', type=int, default=1, help='Starting index to look for galaxy.')
    parser.add_argument('--end', type=int, default=1000, help='Ending index to look for galaxy.')
    parser.add_argument('--maxz',type=float,default=2.0, help='Maximum redshift displayed in plot.')
    parser.add_argument('--NUT', type=bool, default=True, help='If True, it looks for NUT as the most massive galaxy (stars) in the last snapshot.')
    args = parser.parse_args()

    if not isinstance(args.comp, list):
        args.comp = [args.comp]

    # Firstly collect data for control simulation
    if args.model[0] != '/':
        simfolder = os.path.join(os.getcwd(), args.model)
    else:
        simfolder = args.model
        args.model = args.model.split('/')[-1]
        
    if not os.path.exists(simfolder):
            raise Exception('The given control simulation is not found in this directory!')
    
    groupspath = os.path.join(simfolder, 'Groups')

    files = os.listdir(groupspath)

    ozyfiles = []

    for f in files:
        if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
            ozyfiles.append(f)
    ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)

    control_masses = {}
    control_time = []

    for mass_var in args.var:
        control_masses[mass_var] = []
    
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

        for mass_var in args.var:
            m = sim.galaxies[progind].mass[mass_var]
            m = YTQuantity(m, 'code_mass', registry=sim.unit_registry)
            control_masses[mass_var].append(m.in_units('Msun').d)
        control_time.append(thubble)
        
        try:
            progind = sim.galaxies[progind].progen_galaxy_star
        except:
            progind = -1

    # Now, we do the same for the comparison runs

    comp_masses = {}
    comp_times = {}

    for i in range(0, len(args.comp)):

        comp_name = args.comp[i]

        if comp_name[0] != '/':
            simfolder = os.path.join(os.getcwd(), comp_name)
        else:
            simfolder = comp_name
            comp_name = comp_name.split('/')[-1].split('_')[0] + ' '+comp_name.split('/')[-1].split('_')[1]

        comp_masses[comp_name] = {}
        comp_times[comp_name] = []

        if not os.path.exists(simfolder):
            raise Exception('The given simulation %s is not found in this directory!'%comp_name)
    
        groupspath = os.path.join(simfolder, 'Groups')

        files = os.listdir(groupspath)

        ozyfiles = []

        for f in files:
            if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                ozyfiles.append(f)
        ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)

        for mass_var in args.var:
            comp_masses[comp_name][mass_var] = []

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

            for mass_var in args.var:
                m = sim.galaxies[progind].mass[mass_var]
                m = YTQuantity(m, 'code_mass', registry=sim.unit_registry)
                comp_masses[comp_name][mass_var].append(m.in_units('Msun').d)
            comp_times[comp_name].append(thubble)
            
            try:
                progind = sim.galaxies[progind].progen_galaxy_star
            except:
                progind = -1

    # Now, we proceed to do the binning and averaging
    control_time = np.asarray(control_time)
    time_bins = np.arange(control_time.min(),control_time.max(),args.bin/1e+3)

    bin_index = np.digitize(control_time,bins=time_bins,right=True)

    for mass_var in args.var:
        data = np.asarray(control_masses[mass_var])
        control_masses[mass_var] = np.asarray([data[bin_index == i].mean() for i in range(1, len(time_bins))])
    
    for i in range(0, len(args.comp)):

        comp_name = args.comp[i]

        if comp_name[0] != '/':
            simfolder = os.path.join(os.getcwd(), comp_name)
        else:
            simfolder = comp_name
            comp_name = comp_name.split('/')[-1].split('_')[0] + ' '+comp_name.split('/')[-1].split('_')[1]

        t = np.asarray(comp_times[comp_name])
        bin_index = np.digitize(t,bins=time_bins,right=True)

        for mass_var in args.var:
            data = np.asarray(comp_masses[comp_name][mass_var])
            comp_masses[comp_name][mass_var] = np.asarray([data[bin_index == i].mean() for i in range(1, len(time_bins))])
            comp_masses[comp_name][mass_var] = comp_masses[comp_name][mass_var]/control_masses[mass_var]


    # Now add all to a simple plot

    t = 0.5*(time_bins[1:]+time_bins[:-1])

    sub_names = {"stellar":"*","dm":"DM","gas":"gas"}
    marker = itertools.cycle((',', '+', '.', 'o', '*'))


    fig, axes = plt.subplots(len(args.var), 1, sharex=True, figsize=(5,9), dpi=100, facecolor='w', edgecolor='k')

    for i in range(0, len(args.var)):
        ax = axes[i]
        ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
        label = r'$M_{\rm %s}/M_{%s,{\rm %s}}$'%(sub_names[args.var[i]],sub_names[args.var[i]],args.model)
        ax.set_ylabel(label, fontsize=16)
        ax.tick_params(labelsize=12)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")

        for model in args.comp:
            masked_mass= np.ma.masked_invalid(comp_masses[model][args.var[i]])
            ax.plot(t, masked_mass, marker='o', markersize=2, label=model)

        if i == 0:
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
    fig.savefig(os.getcwd()+'/mass_ratio_'+str(args.ind)+'.png', format='png', dpi=200)


