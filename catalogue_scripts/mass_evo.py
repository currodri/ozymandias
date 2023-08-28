"""
This is a simple code that tracks the mass evolution of a particular galaxy across the set of 
snapshots selected, comparing multiple simulations

By: F. Rodriguez Montero (23/02/2021)
Updated: F. Rodriguez Montero (07/02/2022)
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
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "Computer Modern Roman",
})

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
        'cosmoNUTcrmhd\_3e29':'3e29CRMHD',
        'cosmoNUTrticrmhd':'RTCRiMHD'
        }

if __name__ == '__main__':

    # Parse the command line arguments.
    parser = argparse.ArgumentParser(description='Tracking of galaxy gas mass across cosmic time.')
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

    if args.type == 'basic':
        # Now add all to a simple plot
        fig, axes = plt.subplots(len(args.vars),1, sharex=True, figsize=(6,8), dpi=100, facecolor='w', edgecolor='k')
        
        mass_data = {}
        time_data = {}
        for i in range(0,len(args.vars)):
            mass_data[args.vars[i]] = {}
            for j in range(0, len(args.model)):
                if args.model[j][0] != '/':
                    simfolder = os.path.join(os.getcwd(), args.model[i])
                    modelname = args.model[j].replace('_','\_')
                else:
                    simfolder = args.model[j]
                    modelname = args.model[j].split('/')[-1]
                    modelname = args.model[j].replace('_','\_')
                mass_data[args.vars[i]][modelname] = []
                time_data[modelname] = []

        for i in range(0, len(args.model)):
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
                modelname = args.model[i].replace('_','\_')
            else:
                simfolder = args.model[i]
                modelname = args.model[i].split('/')[-1]
                modelname = args.model[i].replace('_','\_')
            print(modelname)
            
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
                m_star = sim.galaxies[progind].mass['stellar']
                for m_var in args.vars:
                    mass_data[m_var][modelname].append(sim.galaxies[progind].mass[m_var].in_units('Msun').d)
                galaxy_masses.append(m_star.in_units('Msun').d)
                galaxy_time.append(thubble)
                bad_mass = False
                if len(galaxy_time) >= 2:
                    if (galaxy_masses[-1]-galaxy_masses[-2])/galaxy_masses[-2] > 0.3:
                        bad_mass = True
                    else:
                        bad_mass = False
                
                if bad_mass == True:
                    print('Deleting')
                    for m_var in args.vars:
                        mass_data[m_var][modelname][-2] = mass_data[m_var][modelname][-1]
                    galaxy_masses[-2] = galaxy_masses[-1]

                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    print('I have lost this galaxy in the snapshot %s'%ozyfile)
                    progind = -1

            time_data[modelname] = galaxy_time
        
        for i in range(0, len(args.vars)):
            ax = axes[i]
            ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
            ax.set_ylabel(r'$M_{\rm %s}/M_{\odot}$'%args.vars[i], fontsize=16)
            ax.set_ylim([8e+6,8e+10])
            ax.set_yscale('log')
            ax.tick_params(labelsize=12)
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.minorticks_on()
            ax.tick_params(which='both',axis="both",direction="in")

            for j in range(0, len(args.model)):
                if args.model[j][0] != '/':
                    simfolder = os.path.join(os.getcwd(), args.model[i])
                    modelname = args.model[j].replace('_','\_')
                else:
                    simfolder = args.model[j]
                    modelname = args.model[j].split('/')[-1]
                    modelname = args.model[j].replace('_','\_')
                ax.plot(time_data[modelname][::-1], mass_data[args.vars[i]][modelname][::-1], marker='o', markersize=2, label=modelname)
        
        axes[0].legend(loc='best', fontsize=14,frameon=False)
        # Add top ticks for redshift
        axR = axes[0].twiny()
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

        fig.subplots_adjust(top=0.93, bottom=0.07,right=0.98,hspace=0.0)
        if args.NUT:
            args.ind = 'NUT'
        fig.savefig(os.getcwd()+'/mass_evolution_'+str(args.ind)+'.png', format='png', dpi=200)

    elif args.type == 'stellar+gas':
        # Now add all to a simple plot
        fig, axes = plt.subplots(1, 2, sharex=True, figsize=(10,5), dpi=100, facecolor='w', edgecolor='k')

        for i in range(0, len(args.model)):
            
            if args.model[i][0] != '/' and args.model[i][0] != '..':
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
            galaxy_gas = []
            galaxy_time = []


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

                m = sim.galaxies[progind].mass['stellar']
                galaxy_masses.append(m.in_units('Msun').d)
                m = sim.galaxies[progind].mass['gas']
                galaxy_gas.append(m.in_units('Msun').d)
                galaxy_time.append(thubble)
                bad_mass = False
                if len(galaxy_time) >= 2:
                    if (galaxy_masses[-1]-galaxy_masses[-2])/galaxy_masses[-2] > 0.3:
                        bad_mass = True
                    else:
                        bad_mass = False
                
                if bad_mass == True:
                    print('Deleting ', ozyfile)
                    galaxy_masses[-2] = galaxy_masses[-1]
                    galaxy_gas[-2] = galaxy_gas[-1]

                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    print('I have lost this galaxy in the snapshot %s'%ozyfile)
                    progind = -1

            stellar_mass = galaxy_masses[::-1]
            gas_mass = galaxy_gas[::-1]
            time = galaxy_time[::-1]
            axes[0].plot(time, stellar_mass, marker='o', markersize=2, label=names[args.model[i]], color=line_dict[args.model[i]])
            axes[1].plot(time, gas_mass, marker='o', markersize=2, label=names[args.model[i]], color=line_dict[args.model[i]])
            
        ax = axes[0]
        ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
        ax.set_ylabel(r'$M_{*}$ [$M_{\odot}$]', fontsize=16)
        ax.set_ylim([9e+7,1e+11])
        ax.set_yscale('log')
        ax.tick_params(labelsize=15)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
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
        axR.tick_params(labelsize=15)

        ax = axes[1]
        ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
        ax.set_ylabel(r'$M_{\rm gas}$ [$M_{\odot}$]', fontsize=16)
        ax.yaxis.set_label_position("right")
        ax.set_ylim([9e7,1e+11])
        ax.set_yscale('log')
        ax.tick_params(labelsize=15,labelleft=False,labelright=True)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
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
        axR.tick_params(labelsize=15)
    
        ax.legend(loc='best', fontsize=14,frameon=False)
        fig.subplots_adjust(top=0.90, bottom=0.1,right=0.93,left=0.08,hspace=0.0,wspace=0.0)
        if args.NUT:
            args.ind = 'NUT'
        fig.savefig(os.getcwd()+'/stellar+gas_evolution_'+str(args.ind)+'.png', format='png', dpi=200)

    elif args.type == 'stellar+gas+cold':
        # Now add all to a simple plot
        fig, axes = plt.subplots(1, 2, sharex=True, figsize=(10,5), dpi=300, facecolor='w', edgecolor='k')

        for i in range(0, len(args.model)):
            
            if args.model[i][0] != '/' and args.model[i][0] != '..':
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
            galaxy_gas = []
            galaxy_cold = []
            galaxy_time = []


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

                m = sim.galaxies[progind].mass['stellar']
                galaxy_masses.append(m.in_units('Msun').d)
                m = sim.galaxies[progind].mass['gas']
                galaxy_gas.append(m.in_units('Msun').d)
                m = sim.galaxies[progind].mass['gas_cold']
                galaxy_cold.append(m.in_units('Msun').d)
                galaxy_time.append(thubble)
                bad_mass = False
                if len(galaxy_time) >= 2:
                    if (galaxy_masses[-1]-galaxy_masses[-2])/galaxy_masses[-2] > 0.3:
                        bad_mass = True
                    else:
                        bad_mass = False
                
                if bad_mass == True:
                    print('Deleting ', ozyfile)
                    galaxy_masses[-2] = galaxy_masses[-1]
                    galaxy_gas[-2] = galaxy_gas[-1]
                    galaxy_cold[-2] = galaxy_cold[-1]

                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    print('I have lost this galaxy in the snapshot %s'%ozyfile)
                    progind = -1

            stellar_mass = galaxy_masses[::-1]
            gas_mass = galaxy_gas[::-1]
            cold_mass = galaxy_cold[::-1]
            time = galaxy_time[::-1]
            axes[0].plot(time, stellar_mass, marker='o', markersize=2, label=names[args.model[i]], color=line_dict[args.model[i]])
            axes[1].plot(time, gas_mass, marker='o', markersize=2, label=names[args.model[i]], color=line_dict[args.model[i]])
            axes[1].plot(time, cold_mass, marker='v', markersize=2, linestyle='--', color=line_dict[args.model[i]], alpha=0.6)
            
        ax = axes[0]
        ax.set_xlabel(r'$t$ [Gyr]', fontsize=20)
        ax.set_ylabel(r'$M_{*}$ [$M_{\odot}$]', fontsize=20)
        ax.set_ylim([9e+7,1e+11])
        ax.set_yscale('log')
        ax.tick_params(labelsize=17)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
        # Add top ticks for redshift
        axR = ax.twiny()
        maxt = cosmo.age(args.maxz).value
        ax.set_xlim(0.0, maxt)
        axR.set_xlim(0.0, maxt)
        topticks1 = np.array([1.5, 2, 3, 4, 6, 10])
        topticks1 = topticks1[topticks1 >= args.maxz]
        topticks2 = cosmo.age(topticks1).value
        axR.set_xticklabels(topticks1)
        axR.set_xticks(topticks2)
        axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
        axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
        axR.set_xlabel(r'$z$', fontsize=20)
        axR.tick_params(labelsize=17)

        ax = axes[1]
        ax.set_xlabel(r'$t$ [Gyr]', fontsize=20)
        ax.set_ylabel(r'$M_{\rm gas}$ [$M_{\odot}$]', fontsize=20)
        ax.yaxis.set_label_position("right")
        ax.set_ylim([9e7,1e+11])
        ax.set_yscale('log')
        ax.tick_params(labelsize=17,labelleft=False,labelright=True)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
        # Add top ticks for redshift
        axR = ax.twiny()
        maxt = cosmo.age(args.maxz).value
        ax.set_xlim(0.0, maxt)
        axR.set_xlim(0.0, maxt)
        topticks1 = np.array([1.5, 2, 3, 4, 6, 10])
        topticks1 = topticks1[topticks1 >= args.maxz]
        topticks2 = cosmo.age(topticks1).value
        axR.set_xticklabels(topticks1)
        axR.set_xticks(topticks2)
        axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
        axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
        axR.set_xlabel(r'$z$', fontsize=20)
        axR.tick_params(labelsize=17)
    
        axes[0].legend(loc='best', fontsize=14,frameon=False, ncol=len(args.model))
        fig.subplots_adjust(top=0.87, bottom=0.13,right=0.91,left=0.08,hspace=0.0,wspace=0.0)
        if args.NUT:
            args.ind = 'NUT'
        fig.savefig(os.getcwd()+'/stellar+gas+cold_evolution_'+str(args.ind)+'.pdf', format='pdf', dpi=300)