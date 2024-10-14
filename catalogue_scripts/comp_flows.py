# Import required libraries
import ozy
from ozy.outflow_inflow import compute_flows, get_flow_name
from ozy.plot_settings import plotting_dictionary
from scipy.signal import medfilt
from scipy.integrate import cumtrapz
import numpy as np
import os
import argparse
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM, z_at_value
import seaborn as sns
from scipy import stats
import matplotlib.patheffects as pe
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
sns.set_theme(style="white")
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "Computer Modern Roman",
})

names = {'cosmoNUThd':'HD',
        'cosmoNUThd\_all\_cr10':'HDcr10',
        'cosmoNUThd\_all\_cr20':'HDcr20',
        'cosmoNUTmhd':'MHD',
        'cosmoNUTcrmhd':'CRMHD',
        'cosmoNUTcrmhd\_nost':'nsCRMHD',
        'cosmoNUTcrmhd\_noheat':'nhCRMHD',
        'cosmoNUTrtcrmhd':'RTCRMHD',
        'cosmoNUTcrmhd\_3e29':'3e29CRMHD'
        }
line_dict = {'cosmoNUThd':'royalblue',
                'cosmoNUThd\_all\_cr10':'lightseagreen',
                'cosmoNUThd\_all\_cr20':'lightcoral',
                'cosmoNUTmhd':'m',
                'cosmoNUTcrmhd':'g',
                'cosmoNUTcrmhd\_nost':'olive',
                'cosmoNUTcrmhd\_noheat':'darkgoldenrod',
                'cosmoNUTcrmhd\_3e29':'firebrick',
                'cosmoNUTrtcrmhd':'firebrick'}

phase_labels = {'hot':'Hot','warm_ionised':'Warm ionised','warm_neutral':'Warm neutral','cold':'Cold'}
phase_colours = {'hot':'maroon','warm_ionised':'orangered','warm_neutral':'orange','cold':'b'}

if __name__ == '__main__':

    # Parse the command line arguments.
    parser = argparse.ArgumentParser(description='Tracking of galaxy energies across cosmic time.')
    parser.add_argument('model', type=str, nargs='+', help='Simulation names from which extract OZY results.')
    parser.add_argument('--type', type=str, default='outflow', help='Flow variable to be plotted.')
    parser.add_argument('--ind', type=int, default=0, help='Index of the galaxy in the ending snapshot.')
    parser.add_argument('--var', type=str, default='massflow_rate', help='Flow variable to be plotted.')
    parser.add_argument('--r', type=float, default=0.2, help='Fraction of rvir at which flow measured.')
    parser.add_argument('--dr', type=float, default=0.01, help='Thickness of the shell in units of rvir.')
    parser.add_argument('--start', type=int, default=10, help='Starting index to look for galaxy.')
    parser.add_argument('--end', type=int, default=1000, help='Ending index to look for galaxy.')
    parser.add_argument('--nbins',type=int,default=1000,help='Number of bins for the underlying PDF.')
    parser.add_argument('--maxz',type=float,default=1.5, help="Maximum redshift displayed in plot.")
    parser.add_argument('--flowtype', type=str, default='outflow', help='Flow type to be plotted.')
    parser.add_argument('--sfr', type=str, default='100Myr',help='What SFR indicator bin to use.')
    parser.add_argument('--NUT', type=bool, default=True, help='If True, it looks for NUT as the most massive galaxy (stars) in the last snapshot.')
    parser.add_argument('--recompute',action='store_true', help='If present, it recomputes 2D profiles.')
    parser.add_argument('--rm_subs',action='store_true', help='If present, it removes the substructures that intersect the spherical shells.')

    args = parser.parse_args()

    if not isinstance(args.model, list):
        args.model = [args.model]
    
    if args.type == 'Kimm+2015':
        # Now add all to a single plot
        fig, axes = plt.subplots(3,3, figsize=(13,13),dpi=100,facecolor='w',edgecolor='k')
        plotting_vars = np.array([[['inflow','massflow_rate',1.0], ['outflow','massflow_rate',0.5], ['outflow','massflow_rate',1.0]],
                            [['inflow','v_sphere_r',1.0], ['outflow','v_sphere_r',0.5], ['outflow','v_sphere_r',1.0]],
                            [['inflow','metallicity',1.0], ['outflow','metallicity',0.5], ['outflow/inflow','density',1.0]]
                            ])

        line_styles = {'cosmoNUThd':':','cosmoNUTmhd':'--','cosmoNUTcrmhd':'-','cosmoNUTcrmhd\_nost':'--','cosmoNUTcrmhd\_noheat':':'}
        flow_style = {'outflow':'-','inflow':'--'}

        data = np.array([[dict(),dict(),dict()],
                        [dict(),dict(),dict()],
                        [dict(),dict(),dict()]
                        ])
        time = dict()
        for i in range(0,3):
            for j in range(0,3):
                for k in range(0, len(args.model)):
                    if i==2 and j==2:
                        data[i,j][args.model[k]] = [[],[]]
                    else:
                        data[i,j][args.model[k]] = []
        for k in range(0, len(args.model)):
            time[args.model[k]] = []

        for i in range(0, len(args.model)):
            if args.model[i] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
            else:
                simfolder = args.model[i]
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            files = os.listdir(groupspath)
            ozyfiles = []

            for f in files:
                if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                    ozyfiles.append(f)
            
            ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
            galaxy_time = []
            galaxy_masses = []

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
                
                if redshift <= args.maxz:
                    try:
                        progind = sim.galaxies[progind].progen_galaxy_star
                    except:
                        progind = -1
                    continue
                h = sim.simulation.hubble_constant
                cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

                # Age of Universe at this redshift
                thubble = cosmo.age(redshift).value

                galaxy_time.append(thubble)

                gal = sim.galaxies[progind]
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
                    try:
                        for k in range(0,3):
                            for j in range(0,3):
                                data[k,j][args.model[i]][-2] = data[k,j][args.model[i]][-1]
                                for t in range(0,2):
                                    data[k,j][args.model[i]][t][-2] = data[k,j][args.model[i]][t][-1]
                    except:
                        pass
                
                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    progind = -1
                
                gf = {}
                gf['inflow_50'] = compute_flows(gal,os.path.join(groupspath, ozyfile),'inflow',rmin=(0.5-0.01,'rvir'),
                                                rmax=(0.5+0.01,'rvir'),save=False,recompute=args.recompute,
                                                remove_subs=args.rm_subs,pdf_bins=args.nbins)
                gf['inflow_100'] = compute_flows(gal,os.path.join(groupspath, ozyfile),'inflow',rmin=(1.0-0.01,'rvir'),
                                                rmax=(1.0+0.01,'rvir'),save=False,recompute=args.recompute,
                                                remove_subs=args.rm_subs,pdf_bins=args.nbins)
                gf['outflow_50'] = compute_flows(gal,os.path.join(groupspath, ozyfile),'outflow',rmin=(0.5-0.01,'rvir'),
                                                rmax=(0.5+0.01,'rvir'),save=False,recompute=args.recompute,
                                                remove_subs=args.rm_subs,pdf_bins=args.nbins)
                gf['outflow_100'] = compute_flows(gal,os.path.join(groupspath, ozyfile),'outflow',rmin=(1.0-0.01,'rvir'),
                                                rmax=(1.0+0.01,'rvir'),save=False,recompute=args.recompute,
                                                remove_subs=args.rm_subs,pdf_bins=args.nbins)
                
                
                factor = 1
                for k in range(0,3):
                    for j in range(0,3):
                        v = plotting_vars[k,j]
                        flow_type = v[0]
                        variable = v[1]
                        radius = float(v[2])
                        d_key = str(int(100*radius))
                        if flow_type == 'inflow' and variable == 'massflow_rate' or flow_type == 'inflow' and variable == 'v_sphere_r':
                            factor = -1
                        else:
                            factor = 1
                        need_to_loop = False
                        if variable == 'density':
                            flow_type = v[0].split('/')
                            radius = [0.5,1.0]
                            need_to_loop = True
                        if need_to_loop:
                            plt_setting = plotting_dictionary['density']
                            for t in range(0,2):
                                d_key = str(int(100*radius[t]))
                                try :
                                    d = gf[flow_type[t]+'_'+d_key].data[variable+'_'+d_key+'rvir_all'].in_units(plt_setting['units'])
                                    data[k,j][args.model[i]][t].append(d)
                                except:
                                    raise KeyError('This variable is not present in this GalacticFlow object: '+str(gf[flow_type[t]+'_'+d_key]))
                        else:
                            plt_setting = plotting_dictionary[variable]
                            try:
                                try:
                                    d = gf[flow_type+'_'+d_key].data[variable+'_'+d_key+'rvir_all'].in_units(plt_setting['units'])
                                except:
                                    d = gf[flow_type+'_'+d_key].data[variable+'_'+d_key+'rvir_all']
                                data[k,j][args.model[i]].append(d*factor)
                            except:
                                raise KeyError('This variable is not present in this GalacticFlow object: '+str(gf[flow_type+'_'+d_key]))
            time[args.model[i]] = galaxy_time
        
        
        for i in range(0,3):
            for j in range(0,3):
                v = plotting_vars[i,j]
                ax = axes[i,j]
                plt_setting = plotting_dictionary[v[1]]
                ax.set_xlabel(r'$t$ [Gyr]', fontsize=12)
                if v[1] != 'v_sphere_r' and v[1] != 'massflow_rate':
                    ax.set_ylabel(plt_setting['label_log'], fontsize=16)
                else:
                    ax.set_ylabel(plt_setting['label'], fontsize=16)
                ax.tick_params(labelsize=12)
                ax.xaxis.set_ticks_position('both')
                ax.yaxis.set_ticks_position('both')
                ax.minorticks_on()
                ax.tick_params(which='both',axis="both",direction="in")
                for k in range(0, len(args.model)):
                    if args.model[k] != '/':
                        temp_name = args.model[k].replace('_','\_')
                    else:
                        temp_name = args.model[k].split('/')[-1]
                        temp_name = temp_name.replace('_','\_')
                    if i==2 and j==2:
                        radius = [0.5,1.0]
                        for t in range(0,2):
                            label = temp_name.split('cosmoNUT')[1] + '('+v[0].split('/')[t]+','+str(radius[t])+r'$r_{\rm vir,DM}$)'
                            median = medfilt(data[i,j][args.model[k]][t],7)
                            ax.plot(time[args.model[k]], np.log10(median),label=label,linestyle=flow_style[v[0].split('/')[t]],color=line_dict[temp_name])
                            ax.set_ylim([-28,-22])
                    else:
                        if v[1] != 'v_sphere_r' and v[1] != 'massflow_rate':
                            #median = RunningMedian(data[i,j][args.model[k]],3)
                            median = medfilt(data[i,j][args.model[k]],5)
                            ax.plot(time[args.model[k]], np.log10(median),linestyle=line_styles[temp_name],color=line_dict[temp_name])
                        else:
                            median = medfilt(data[i,j][args.model[k]],5)
                            ax.plot(time[args.model[k]], median,label=temp_name,linestyle=line_styles[temp_name],color=line_dict[temp_name])
                        

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
                axR.set_xlabel(r'$z$', fontsize=10)
                axR.tick_params(labelsize=10)
                if i==0 and j==0:
                    ax.legend(loc='lower right', fontsize=14,frameon=False)
                if i==2 and j==2:
                    ax.legend(loc='upper left', fontsize=14,frameon=False)
        for i in range(0,3):
            axes[i,0].text(0.05, 0.95, r'inflows at $r=R_{\rm vir,DM}$',
                        transform=axes[i,0].transAxes, fontsize=17,verticalalignment='top',
                        color='black')
        for i in range(0,3):
            axes[i,1].text(0.05, 0.95, r'outflows at $r=0.5R_{\rm vir,DM}$',
                        transform=axes[i,1].transAxes, fontsize=17,verticalalignment='top',
                        color='black')
        for i in range(0,2):
            axes[i,2].text(0.05, 0.95, r'outflows at $r=R_{\rm vir,DM}$',
                        transform=axes[i,2].transAxes, fontsize=17,verticalalignment='top',
                        color='black')
        ax = axes[2,2]
        time = np.linspace(np.min(time[args.model[0]]),np.max(time[args.model[0]]),100)
        redshift = [z_at_value(cosmo.age,t*u.Gyr) for t in time]
        rho_crit = cosmo.critical_density(redshift)
        rho_halo = 178*rho_crit.value*sim.simulation.omega_baryon/sim.simulation.omega_matter
        ax.plot(time,np.log10(rho_halo),linewidth=4,color='k',alpha=0.5,label=r'$\rho_{\rm b}$')
        ax.legend(loc='best', fontsize=14,frameon=False)

        if args.NUT:
            args.ind = 'NUT'
        fig.subplots_adjust(hspace=0.2,wspace=0.2,left=0.065,right=0.98,top=0.95,bottom=0.05)
        fig.savefig(os.getcwd()+'/flowcomp_Kimm+2015_'+args.var+'_'+str(args.ind)+'.png', format='png', dpi=200)
    
    elif args.type == 'cumulative':
        # Now add all to a single plot
        fig, ax = plt.subplots(1,1, figsize=(7,5),dpi=100,facecolor='w',edgecolor='k')
        plt_setting = plotting_dictionary[args.var]
        ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
        ax.set_ylabel(r'$\log{M/M_{\odot}}$', fontsize=16)
        ax.tick_params(labelsize=12)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
        ax.set_yscale('log')
        for i in range(0, len(args.model)):
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            files = os.listdir(groupspath)
            ozyfiles = []

            for f in files:
                if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                    ozyfiles.append(f)
            
            ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
            galaxy_time = []
            galaxy_masses = []
            flow_values = []

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
                if redshift <= args.maxz:
                    try:
                        progind = sim.galaxies[progind].progen_galaxy_star
                    except:
                        progind = -1
                    continue
                h = sim.simulation.hubble_constant
                cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

                # Age of Universe at this redshift
                thubble = cosmo.age(redshift).value

                galaxy_time.append(thubble)

                gal = sim.galaxies[progind]
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
                    try:
                        flow_values[-2] = flow_values[-1]
                    except:
                        pass
                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    progind = -1
                
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.5*args.dr,'rvir'),
                                    rmax=(args.r+0.5*args.dr,'rvir'),save=True,recompute=args.recompute,
                                    remove_subs=args.rm_subs,pdf_bins=args.nbins)
                
                d_key = str(int(100*args.r))
                try:
                    d = gf.data[args.var+'_'+d_key+'rvir_all'].in_units(plt_setting['units'])
                    if args.flowtype == 'inflow':
                        d = -d
                    flow_values.append(d)
                except:
                    raise KeyError('This variable is not present in this GalacticFlow object: '+str(gf.data[args.var+'_'+d_key+'rvir']))

            # Integrate curve
            print(flow_values)
            flow_values = cumtrapz(np.asarray(flow_values), np.asarray(galaxy_time)*1e+9)
            print(flow_values)
            galaxy_time = np.asarray(galaxy_time)
            galaxy_time = 0.5*(galaxy_time[:-1]+galaxy_time[1:])
            ax.plot(galaxy_time,-1*flow_values[::-1],marker='o', markersize=2, label=args.model[i])

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
        ax.legend(loc='best', fontsize=14,frameon=False)

        if args.NUT:
            args.ind = 'NUT'
        fig.savefig(os.getcwd()+'/cumulative_'+args.flowtype+'_'+args.var+'_'+str(args.ind)+'_'+d_key+'rvir.png', format='png', dpi=200)
    
    elif args.type == 'separate_phases':
        fig, axes = plt.subplots(3,len(args.model), figsize=(13,13),dpi=100,facecolor='w',edgecolor='k', sharex='col', sharey='row')

        if args.flowtype == 'outflow':
            line_styles = {'hot':'-','warm_ionised':'--','warm_neutral':'-.','cold':':'}
            labels = {'hot':r'$T > 10^5$K','warm_ionised':r'$9\cdot 10^3 \text{K} < T < 10^5$K','warm_neutral':r'$1\cdot 10^3 \text{K} < T < 9\cdot 10^3$K','cold':r'$1\cdot 10^3 \text{K} < T$'}
        elif args.flowtype == 'inflow':
            line_styles = {'hot':'-','warm_ionised':'--','warm_neutral':'-.','cold':':'}
        plotting_vars = np.array([['massflow_rate',args.r],['v_sphere_r',args.r],['metallicity',args.r]])

        for i in range(0,len(args.model)):
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            files = os.listdir(groupspath)
            ozyfiles = []

            for f in files:
                if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                    ozyfiles.append(f)
            
            ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
            galaxy_time = []
            galaxy_masses = []

            progind = args.ind
            if args.NUT:
                sim = ozy.load(os.path.join(groupspath, ozyfiles[0]))
                virial_mass = [i.virial_quantities['mass'] for i in sim.galaxies]
                progind = np.argmax(virial_mass)

            data = np.array([dict(),dict(),dict()])
            for d in range(0,3):
                for f in line_styles.keys():
                    data[d][f] = []
            galaxy_time = []
            for ozyfile in ozyfiles:
                if progind == -1:
                    continue
                # Load OZY file
                sim = ozy.load(os.path.join(groupspath, ozyfile))
                # Initialise simulation parameters
                redshift = sim.simulation.redshift
                
                if redshift <= args.maxz:
                    try:
                        progind = sim.galaxies[progind].progen_galaxy_star
                    except:
                        progind = -1
                    continue
                h = sim.simulation.hubble_constant
                cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

                # Age of Universe at this redshift
                thubble = cosmo.age(redshift).value

                galaxy_time.append(thubble)

                gal = sim.galaxies[progind]
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
                    try:
                        for d in range(0,3):
                            for f in line_styles:
                                data[d][f][-2] = data[d][f][-1]
                    except:
                        pass
                
                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    progind = -1
                
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.5*args.dr,'rvir'),
                                    rmax=(args.r+0.5*args.dr,'rvir'),save=True,recompute=args.recompute,separate_phases=True,
                                    remove_subs=args.rm_subs,pdf_bins=args.nbins)
                
                d_key = str(int(100*args.r))
                factor = 1
                for d in range(0,3):
                        for f in line_styles.keys():
                            v = plotting_vars[d][0]
                            plt_setting = plotting_dictionary[v]
                            print(v+'_'+d_key+'rvir_'+str(f))
                            try:
                                value = gf.data[v+'_'+d_key+'rvir_'+str(f)].in_units(plt_setting['units'])
                                if args.flowtype == 'inflow':
                                    value = -value
                                data[d][f].append(value)
                            except:
                                raise KeyError('This variable is not present in this GalacticFlow object: '+str(gf.data[v+'_'+d_key+'rvir_'+str(f)]))

            for d in range(0,3):
                v = plotting_vars[d][0]
                ax = axes[d,i]
                plt_setting = plotting_dictionary[v]
                ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
                if i == 0:
                    if v != 'v_sphere_r' and v != 'massflow_rate':
                        ax.set_ylabel(plt_setting['label_log'], fontsize=16)
                    else:
                        ax.set_ylabel(plt_setting['label'], fontsize=16)
                ax.tick_params(labelsize=12)
                ax.xaxis.set_ticks_position('both')
                ax.yaxis.set_ticks_position('both')
                ax.minorticks_on()
                ax.tick_params(which='both',axis="both",direction="in")
                for f in line_styles.keys():
                    label = labels[f]
                    median = medfilt(data[d][f],7)
                    if v != 'v_sphere_r' and v != 'massflow_rate':
                        ax.plot(galaxy_time, np.log10(median),linestyle=line_styles[f],color=line_dict[args.model[i]],label=label)
                    else:
                        ax.plot(galaxy_time, median,linestyle=line_styles[f],color=line_dict[args.model[i]],label=label)
                if d == 0:
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
                    axR.tick_params(labelsize=14)

                    ax.legend(loc='best', fontsize=14,frameon=False)
        if args.NUT:
            args.ind = 'NUT'
        fig.subplots_adjust(hspace=0.0,wspace=0.0)#,left=0.065,right=0.98,top=0.95,bottom=0.05)
        fig.savefig(os.getcwd()+'/'+args.flowtype+'_phasecomp_'+args.var+'_'+str(args.ind)+'_'+d_key+'rvir.png', format='png', dpi=200)
    
    elif args.type == 'instantaneous_phases':
        nrow = int(1 + len(args.model))
        height_ratios = []
        height_ratios.append(1)
        for i in range(0, len(args.model)):
            height_ratios.append(0.3)
        figsize = plt.figaspect(float((5.0 * len(args.model)) / (5.0 * 2)))
        fig = plt.figure(figsize=2*figsize, facecolor='w', edgecolor='k')
        plot_grid = fig.add_gridspec(nrow, 1, wspace=0, hspace=0, height_ratios=height_ratios)
        axes = []
        axes.append(fig.add_subplot(plot_grid[0]))
        for i in range(1,nrow):
            axes.append(fig.add_subplot(plot_grid[i],sharex=axes[0]))
        axes = np.asarray(axes)

        # Get data
        plt_setting = plotting_dictionary[args.var]
        for i in range(0, len(args.model)):
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
                args.model[i] = args.model[i].replace('_','\_')
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
                args.model[i] = args.model[i].replace('_','\_')
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            files = os.listdir(groupspath)
            ozyfiles = []

            for f in files:
                if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                    ozyfiles.append(f)
            
            ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
            galaxy_time = []
            galaxy_masses = []
            flow_values = {'all':[],'hot':[],'warm_ionised':[],'warm_neutral':[],'cold':[]}

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
                if redshift <= args.maxz:
                    try:
                        progind = sim.galaxies[progind].progen_galaxy_star
                    except:
                        progind = -1
                    continue
                h = sim.simulation.hubble_constant
                cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

                # Age of Universe at this redshift
                thubble = cosmo.age(redshift).value

                galaxy_time.append(thubble)

                gal = sim.galaxies[progind]
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
                    try:
                        flow_values[-2] = flow_values[-1]
                    except:
                        pass
                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    progind = -1
                
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.5*args.dr,'rvir'),
                                    rmax=(args.r+0.5*args.dr,'rvir'),save=True,recompute=args.recompute,
                                    remove_subs=args.rm_subs,pdf_bins=args.nbins)
                
                d_key = str(int(100*args.r))
                try:
                    for ftype in flow_values.keys():
                        d = gf.data[args.var+'_'+d_key+'rvir_'+ftype].in_units(plt_setting['units'])
                        if args.flowtype == 'inflow':
                            d = -d
                        flow_values[ftype].append(d)
                except:
                    raise KeyError('This variable is not present in this GalacticFlow object: '+str(gf.data[args.var+'_'+d_key+'rvir']))
            
            # Convert to numpy arrays
            for k in flow_values.keys():
                flow_values[k] = np.asarray(flow_values[k])

            # Add total to main axes
            axes[0].step(galaxy_time[::-1],flow_values['all'][::-1], where='mid',alpha=0.4,color=line_dict[args.model[i]])
            all_median = medfilt(flow_values['all'][::-1],5)
            axes[0].plot(galaxy_time[::-1],all_median,marker='o', markersize=4, label=args.model[i],color=line_dict[args.model[i]])

            # Now compute phase contributions and add to specific plots
            hot = 100*(flow_values['hot'][::-1]/flow_values['all'][::-1])
            warm_ionised = 100*((flow_values['hot'][::-1]+flow_values['warm_ionised'][::-1])/flow_values['all'][::-1])
            warm_neutral = 100*((flow_values['hot'][::-1]+flow_values['warm_ionised'][::-1]+flow_values['warm_neutral'][::-1])/flow_values['all'][::-1])
            cold = 100*((flow_values['hot'][::-1]+flow_values['warm_ionised'][::-1]+flow_values['warm_neutral'][::-1]+flow_values['cold'][::-1])/flow_values['all'][::-1])
            
            print(args.model[i], cold[0],warm_neutral[0],warm_ionised[0],hot[0])
            axes[1+i].plot(galaxy_time[::-1],cold,marker='o', markersize=2, label=phase_labels['cold'],color=phase_colours['cold'])
            axes[1+i].fill_between(galaxy_time[::-1],cold,color=phase_colours['cold'],alpha=0.2)
            
            axes[1+i].plot(galaxy_time[::-1],warm_neutral,marker='o', markersize=2, label=phase_labels['warm_neutral'],color=phase_colours['warm_neutral'])
            axes[1+i].fill_between(galaxy_time[::-1],warm_neutral,color=phase_colours['warm_neutral'],alpha=0.2)
            
            axes[1+i].plot(galaxy_time[::-1],warm_ionised,marker='o', markersize=2, label=phase_labels['warm_ionised'],color=phase_colours['warm_ionised'])
            axes[1+i].fill_between(galaxy_time[::-1],warm_ionised,color=phase_colours['warm_ionised'],alpha=0.2)

            axes[1+i].plot(galaxy_time[::-1],hot,marker='o', markersize=2, label=phase_labels['hot'],color=phase_colours['hot'])
            axes[1+i].fill_between(galaxy_time[::-1],hot,color=phase_colours['hot'],alpha=0.2)
            axes[1+i].text(0.65, 0.3, args.model[i],
                        transform=axes[1+i].transAxes, fontsize=16,verticalalignment='top',
                        color='black')
            axes[1+i].tick_params(labelsize=14)
            axes[1+i].xaxis.set_ticks_position('both')
            axes[1+i].yaxis.set_ticks_position('both')
            axes[1+i].minorticks_on()
            axes[1+i].tick_params(which='both',axis="both",direction="in",bottom=False)

        # Add legends to axes
        ax = axes[0]
        ax.set_ylabel(r'$\dot{M}_{\rm %s}$ [$M_{\odot}$/yr]'%args.flowtype, fontsize=16)
        ax.tick_params(labelsize=14)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in",bottom=False)
        ax.set_yscale('log')
        ax.legend(loc='best',fontsize=16,frameon=False)
        ax.text(0.05, 0.95, r'$r= %.1f R_{\rm vir,DM}$'%(float(args.r)),
                        transform=ax.transAxes, fontsize=16,verticalalignment='top',
                        color='black')

        # Add top ticks for redshift
        axR = ax.twiny()
        maxt = cosmo.age(args.maxz).value
        ax.set_xlim(0.0, maxt)
        axR.set_xlim(0.0, maxt)
        topticks1 = np.array([1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 13.0])
        topticks1 = topticks1[topticks1 >= args.maxz]
        topticks2 = cosmo.age(topticks1).value
        axR.set_xticklabels(topticks1)
        axR.set_xticks(topticks2)
        axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
        axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
        axR.set_xlabel(r'$z$', fontsize=18)
        axR.tick_params(labelsize=14)

        ax = axes[1]
        ax.legend(loc='best',fontsize=16,frameon=False,ncol=2)

        # Add label for rest of axes
        ax = axes[-1]
        ax.set_xlabel(r'$t$ [Gyr]', fontsize=18)
        ax.tick_params(which='both',axis="both",direction="in",bottom=True)
        fig.text(0.01, 0.25, r'Mass fraction [$\%$]', fontsize=18,
                rotation='vertical',color='black',va='center')

        if args.NUT:
            args.ind = 'NUT'
        fig.subplots_adjust(top=0.95,bottom=0.05,left=0.1,right=0.98)
        fig.savefig(os.getcwd()+'/instantaneous_'+args.flowtype+'_'+args.var+'_'+str(args.ind)+'_'+d_key+'rvir.png', format='png', dpi=200)
    elif args.type == 'instantaneous_phases_sfr':
        nrow = int(2 + len(args.model))
        height_ratios = []
        height_ratios.append(0.5)
        height_ratios.append(0.5)
        for i in range(0, len(args.model)):
            height_ratios.append(0.3)
        figsize = plt.figaspect(float((5.0 * len(args.model)) / (5.0 * 2)))
        fig = plt.figure(figsize=2*figsize, facecolor='w', edgecolor='k')
        plot_grid = fig.add_gridspec(nrow, 1, wspace=0, hspace=0, height_ratios=height_ratios)
        axes = []
        axes.append(fig.add_subplot(plot_grid[0]))
        for i in range(1,nrow):
            axes.append(fig.add_subplot(plot_grid[i],sharex=axes[0]))
        axes = np.asarray(axes)
        
        phase_labels = {'hot':r'\textbf{Hot}','warm_ionised':r'\textbf{Warm ionised}','warm_neutral':r'\textbf{Warm neutral}','cold':r'\textbf{Cold}'}
        phase_colours = {'hot':'r','warm_ionised':'orange','warm_neutral':'gold','cold':'b'}
        line_dict = {'cosmoNUThd':'royalblue','cosmoNUTmhd':'m','cosmoNUTcrmhd':'g',
                    'cosmoNUTcrmhd\_nost':'olive','cosmoNUTcrmhd\_noheat':'darkgoldenrod',
                    'cosmoNUTcrmhd\_3e29':'olivedrab'}
        # Get data
        plt_setting = plotting_dictionary[args.var]
        for i in range(0, len(args.model)):
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
                args.model[i] = args.model[i].replace('_','\_')
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
                args.model[i] = args.model[i].replace('_','\_')
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            files = os.listdir(groupspath)
            ozyfiles = []

            for f in files:
                if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                    ozyfiles.append(f)
            
            ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
            galaxy_time = []
            galaxy_masses = []
            flow_values = {'all':[],'hot':[],'warm_ionised':[],'warm_neutral':[],'cold':[]}
            sfr = []

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
                print(args.model[i],redshift,ozyfile)
                if redshift <= args.maxz:
                    try:
                        progind = sim.galaxies[progind].progen_galaxy_star
                    except:
                        progind = -1
                    continue
                h = sim.simulation.hubble_constant
                cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

                # Age of Universe at this redshift
                thubble = cosmo.age(redshift).value

                galaxy_time.append(thubble)

                gal = sim.galaxies[progind]
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
                    try:
                        flow_values[-2] = flow_values[-1]
                    except:
                        pass
                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    progind = -1
                
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.5*args.dr,'rvir'),
                                    rmax=(args.r+0.5*args.dr,'rvir'),save=True,recompute=args.recompute,
                                    remove_subs=args.rm_subs,pdf_bins=args.nbins)
                
                d_key = str(int(100*args.r))
                try:
                    for ftype in flow_values.keys():
                        d = gf.data[args.var+'_'+d_key+'rvir_'+ftype].in_units(plt_setting['units'])
                        if args.flowtype == 'inflow':
                            d = -d
                        flow_values[ftype].append(d)
                except:
                    raise KeyError('This variable is not present in this GalacticFlow object: '+str(gf.data[args.var+'_'+d_key+'rvir']))
                
                # And now get the SFR
                print(args.model[i],flow_values['hot'][-1],flow_values['warm_ionised'][-1],flow_values['warm_neutral'][-1],flow_values['cold'][-1])
                sfr.append(gal.sfr[args.sfr].in_units('Msun/yr'))
            
            # Convert to numpy arrays
            for k in flow_values.keys():
                flow_values[k] = np.asarray(flow_values[k])

            # Add total to first axes
            axes[0].step(galaxy_time[::-1],flow_values['all'][::-1], where='mid',alpha=0.4,color=line_dict[args.model[i]])
            all_median = medfilt(flow_values['all'][::-1],5)
            axes[0].plot(galaxy_time[::-1],all_median, linewidth=2.5,label=names[args.model[i]],color=line_dict[args.model[i]])

            # Add SFR to second axes
            axes[1].step(galaxy_time[::-1],sfr[::-1], where='mid',alpha=0.4,color=line_dict[args.model[i]])
            axes[1].plot(galaxy_time[::-1],medfilt(sfr[::-1],5),linewidth=2.5, label=names[args.model[i]],color=line_dict[args.model[i]])

            # Now compute phase contributions and add to specific plots
            hot = 100*(flow_values['hot'][::-1]/flow_values['all'][::-1])
            warm_ionised = 100*((flow_values['hot'][::-1]+flow_values['warm_ionised'][::-1])/flow_values['all'][::-1])
            warm_neutral = 100*((flow_values['hot'][::-1]+flow_values['warm_ionised'][::-1]+flow_values['warm_neutral'][::-1])/flow_values['all'][::-1])
            cold = 100*((flow_values['hot'][::-1]+flow_values['warm_ionised'][::-1]+flow_values['warm_neutral'][::-1]+flow_values['cold'][::-1])/flow_values['all'][::-1])
            
            hot  = medfilt(hot,5)
            warm_ionised = medfilt(warm_ionised,5)
            warm_neutral = medfilt(warm_neutral,5)
            cold = medfilt(cold,5)
            
            axes[2+i].plot(galaxy_time[::-1],cold, label=phase_labels['cold'],color=phase_colours['cold'],linewidth=1.5)
            axes[2+i].fill_between(galaxy_time[::-1],cold,color=phase_colours['cold'],alpha=0.2)
            
            axes[2+i].plot(galaxy_time[::-1],warm_neutral, label=phase_labels['warm_neutral'],color=phase_colours['warm_neutral'],linewidth=1.5)
            axes[2+i].fill_between(galaxy_time[::-1],warm_neutral,color=phase_colours['warm_neutral'],alpha=0.2)
            
            axes[2+i].plot(galaxy_time[::-1],warm_ionised, label=phase_labels['warm_ionised'],color=phase_colours['warm_ionised'],linewidth=1.5)
            axes[2+i].fill_between(galaxy_time[::-1],warm_ionised,color=phase_colours['warm_ionised'],alpha=0.2)

            axes[2+i].plot(galaxy_time[::-1],hot, label=phase_labels['hot'],color=phase_colours['hot'],linewidth=1.5)
            axes[2+i].fill_between(galaxy_time[::-1],hot,color=phase_colours['hot'],alpha=0.2)
            axes[2+i].text(0.77, 0.3, names[args.model[i]],
                        transform=axes[2+i].transAxes, fontsize=16,verticalalignment='top',
                        color='black')
            axes[2+i].tick_params(labelsize=14)
            axes[2+i].xaxis.set_ticks_position('both')
            axes[2+i].yaxis.set_ticks_position('both')
            axes[2+i].minorticks_on()
            axes[2+i].tick_params(which='both',axis="both",direction="in",labelbottom='off')

        # Add legends to axes
        ax = axes[0]
        ax.set_ylim([2e-2,50])
        ax.set_ylabel(r'$\dot{M}_{\rm %s}$ [M$_{\odot}$/yr]'%args.flowtype, fontsize=20)
        ax.tick_params(labelsize=16)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
        ax.set_yscale('log')
        ax.legend(loc='best',fontsize=16,frameon=False)
        ax.text(0.77, 0.95, r'$r= %.1f R_{\rm vir,DM}$'%(float(args.r)),
                        transform=ax.transAxes, fontsize=16,verticalalignment='top',
                        color='black')
        ax = axes[1]
        ax.set_ylabel(r'SFR [M$_{\odot}$/yr]', fontsize=20)
        ax.tick_params(labelsize=16)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
        ax.set_yscale('log')

        # Add top ticks for redshift
        ax = axes[0]
        axR = ax.twiny()
        maxt = cosmo.age(args.maxz).value
        mint = cosmo.age(13.0).value
        ax.set_xlim(mint, maxt)
        axR.set_xlim(mint, maxt)
        topticks1 = np.array([1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
        topticks1 = topticks1[topticks1 >= args.maxz]
        topticks2 = cosmo.age(topticks1).value
        axR.set_xticklabels(topticks1)
        axR.set_xticks(topticks2)
        axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
        axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
        axR.set_xlabel(r'$z$', fontsize=20)
        axR.tick_params(labelsize=16)

        ax = axes[2]
        ax.legend(loc='upper left',fontsize=16,frameon=False,ncol=2)

        # Add label for rest of axes
        ax = axes[-1]
        ax.set_xlabel(r'$t$ [Gyr]', fontsize=20)
        ax.tick_params(which='both',axis="both",direction="in",bottom=True)
        ax.tick_params(labelsize=16)
        fig.text(0.01, 0.25, r'Mass fraction [$\%$]', fontsize=20,
                rotation='vertical',color='black',va='center')

        if args.NUT:
            args.ind = 'NUT'
        fig.subplots_adjust(top=0.95,bottom=0.05,left=0.13,right=0.97)
        if args.rm_subs:
            fig.savefig(os.getcwd()+'/instantaneous_sfr_'+args.flowtype+'_'+args.var+'_'+str(args.ind)+'_'+d_key+'rvir_rmsubs.pdf', format='pdf', dpi=300)
        else:
            fig.savefig(os.getcwd()+'/instantaneous_sfr_'+args.flowtype+'_'+args.var+'_'+str(args.ind)+'_'+d_key+'rvir.pdf', format='pdf', dpi=300)
    elif args.type == 'instantaneous_sfr_eta':
        nrow = 3
        fig = plt.figure(figsize=(6,9), facecolor='w', edgecolor='k')
        plot_grid = fig.add_gridspec(nrow, 1, wspace=0, hspace=0)
        axes = []
        axes.append(fig.add_subplot(plot_grid[0]))
        for i in range(1,nrow):
            axes.append(fig.add_subplot(plot_grid[i],sharex=axes[0]))
        axes = np.asarray(axes)
        
        phase_labels = {'hot':r'\textbf{Hot}','warm_ionised':r'\textbf{Warm ionised}','warm_neutral':r'\textbf{Warm neutral}','cold':r'\textbf{Cold}'}
        phase_colours = {'hot':'r','warm_ionised':'orange','warm_neutral':'gold','cold':'b'}

        # Get data
        plt_setting = plotting_dictionary[args.var]
        for i in range(0, len(args.model)):
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
                args.model[i] = args.model[i].replace('_','\_')
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
                args.model[i] = args.model[i].replace('_','\_')
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            files = os.listdir(groupspath)
            ozyfiles = []

            for f in files:
                if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                    ozyfiles.append(f)
            
            ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
            galaxy_time = []
            galaxy_masses = []
            flow_values = {'all':[],'hot':[],'warm_ionised':[],'warm_neutral':[],'cold':[]}
            sfr = []

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
                print(args.model[i],redshift,ozyfile)
                if redshift <= args.maxz:
                    try:
                        progind = sim.galaxies[progind].progen_galaxy_star
                    except:
                        progind = -1
                    continue
                h = sim.simulation.hubble_constant
                cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

                # Age of Universe at this redshift
                thubble = cosmo.age(redshift).value

                galaxy_time.append(thubble)

                gal = sim.galaxies[progind]
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
                    try:
                        flow_values[-2] = flow_values[-1]
                    except:
                        pass
                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    progind = -1
                
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.5*args.dr,'rvir'),
                                    rmax=(args.r+0.5*args.dr,'rvir'),save=True,recompute=args.recompute,
                                    remove_subs=args.rm_subs,pdf_bins=args.nbins)
                
                d_key = str(int(100*args.r))
                try:
                    for ftype in flow_values.keys():
                        d = gf.data[args.var+'_'+d_key+'rvir_'+ftype].in_units(plt_setting['units'])
                        if args.flowtype == 'inflow':
                            d = -d
                        flow_values[ftype].append(d)
                except:
                    raise KeyError('This variable is not present in this GalacticFlow object: '+str(gf.data[args.var+'_'+d_key+'rvir']))
                
                # And now get the SFR
                print(args.model[i],flow_values['hot'][-1],flow_values['warm_ionised'][-1],flow_values['warm_neutral'][-1],flow_values['cold'][-1])
                sfr.append(gal.sfr[args.sfr].in_units('Msun/yr'))
            
            # Convert to numpy arrays
            for k in flow_values.keys():
                flow_values[k] = np.asarray(flow_values[k])

            # Add total to first axes
            axes[0].step(galaxy_time[::-1],flow_values['all'][::-1], where='mid',alpha=0.4,color=line_dict[args.model[i]])
            all_median = medfilt(flow_values['all'][::-1],5)
            axes[0].plot(galaxy_time[::-1],all_median, linewidth=2.5,label=names[args.model[i]],color=line_dict[args.model[i]])

            # Add SFR to second axes
            axes[1].step(galaxy_time[::-1],sfr[::-1], where='mid',alpha=0.4,color=line_dict[args.model[i]])
            axes[1].plot(galaxy_time[::-1],medfilt(sfr[::-1],5),linewidth=2.5, label=names[args.model[i]],color=line_dict[args.model[i]])

            # Add mass-loading (eta) to third axis
            axes[2].plot(galaxy_time[::-1],all_median/medfilt(sfr[::-1],5),linewidth=2.5, label=names[args.model[i]],color=line_dict[args.model[i]])

        # Add legends to axes
        ax = axes[0]
        ax.set_ylim([2e-2,50])
        ax.set_ylabel(r'$\dot{M}_{\rm %s}$ [M$_{\odot}$/yr]'%args.flowtype, fontsize=20)
        ax.tick_params(labelsize=16)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
        ax.set_yscale('log')
        ax.legend(loc='best',fontsize=16,frameon=False)
        ax.text(0.7, 0.95, r'$r= %.1f R_{\rm vir,DM}$'%(float(args.r)),
                        transform=ax.transAxes, fontsize=16,verticalalignment='top',
                        color='black')
        ax = axes[1]
        ax.set_ylabel(r'SFR [M$_{\odot}$/yr]', fontsize=20)
        ax.tick_params(labelsize=16)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
        ax.set_yscale('log')

        ax = axes[2]
        ax.set_ylabel(r'$\eta$', fontsize=20)
        ax.tick_params(labelsize=16)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
        ax.set_yscale('log')

        # Add top ticks for redshift
        ax = axes[0]
        axR = ax.twiny()
        maxt = cosmo.age(args.maxz).value
        mint = cosmo.age(13.0).value
        ax.set_xlim(mint, maxt)
        axR.set_xlim(mint, maxt)
        topticks1 = np.array([1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
        topticks1 = topticks1[topticks1 >= args.maxz]
        topticks2 = cosmo.age(topticks1).value
        axR.set_xticklabels(topticks1)
        axR.set_xticks(topticks2)
        axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
        axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
        axR.set_xlabel(r'$z$', fontsize=20)
        axR.tick_params(labelsize=16)

        # Add label for rest of axes
        ax = axes[-1]
        ax.set_xlabel(r'$t$ [Gyr]', fontsize=20)
        ax.tick_params(which='both',axis="both",direction="in",bottom=True)
        ax.tick_params(labelsize=16)

        if args.NUT:
            args.ind = 'NUT'
        fig.subplots_adjust(top=0.94,bottom=0.06,left=0.14,right=0.98)
        if args.rm_subs:
            fig.savefig(os.getcwd()+'/instantaneous_sfr_eta_'+args.flowtype+'_'+args.var+'_'+str(args.ind)+'_'+d_key+'rvir_rmsubs.pdf', format='pdf', dpi=300)
        else:
            fig.savefig(os.getcwd()+'/instantaneous_sfr_eta_'+args.flowtype+'_'+args.var+'_'+str(args.ind)+'_'+d_key+'rvir.pdf', format='pdf', dpi=300)
    elif args.type == 'instantaneous_phases_ratio':
        nrow = int(2 + len(args.model))
        height_ratios = []
        height_ratios.append(0.5)
        height_ratios.append(0.5)
        for i in range(0, len(args.model)):
            height_ratios.append(0.3)
        figsize = plt.figaspect(float((5.0 * len(args.model)) / (5.0 * 2)))
        fig = plt.figure(figsize=2*figsize, facecolor='w', edgecolor='k')
        plot_grid = fig.add_gridspec(nrow, 1, wspace=0, hspace=0, height_ratios=height_ratios)
        axes = []
        axes.append(fig.add_subplot(plot_grid[0]))
        for i in range(1,nrow):
            axes.append(fig.add_subplot(plot_grid[i],sharex=axes[0]))
        axes = np.asarray(axes)

        flow_arrays = {}

        maxt = 4.0 # Gyr
        mint = 0.5 # Gyr
        bins = np.linspace(mint,maxt,30)

        # Get data
        plt_setting = plotting_dictionary[args.var]
        for i in range(0, len(args.model)):
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
                args.model[i] = args.model[i].replace('_','\_')
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
                args.model[i] = args.model[i].replace('_','\_')
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            files = os.listdir(groupspath)
            ozyfiles = []

            for f in files:
                if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                    ozyfiles.append(f)
            
            ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
            galaxy_time = []
            galaxy_masses = []
            flow_values = {'all':[],'hot':[],'warm_ionised':[],'warm_neutral':[],'cold':[]}
            
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
                if redshift <= args.maxz:
                    try:
                        progind = sim.galaxies[progind].progen_galaxy_star
                    except:
                        progind = -1
                    continue
                h = sim.simulation.hubble_constant
                cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

                # Age of Universe at this redshift
                thubble = cosmo.age(redshift).value

                galaxy_time.append(thubble)

                gal = sim.galaxies[progind]
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
                    try:
                        flow_values[-2] = flow_values[-1]
                    except:
                        pass
                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    progind = -1
                
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.5*args.dr,'rvir'),
                                    rmax=(args.r+0.5*args.dr,'rvir'),save=True,recompute=args.recompute,
                                    remove_subs=args.rm_subs,pdf_bins=args.nbins)
                
                d_key = str(int(100*args.r))
                try:
                    for ftype in flow_values.keys():
                        d = gf.data[args.var+'_'+d_key+'rvir_'+ftype].in_units(plt_setting['units'])
                        if args.flowtype == 'inflow':
                            d = -d
                        flow_values[ftype].append(d)
                except:
                    raise KeyError('This variable is not present in this GalacticFlow object: '+str(gf.data[args.var+'_'+d_key+'rvir']))
            
            # Convert to numpy arrays
            for k in flow_values.keys():
                flow_values[k] = np.asarray(flow_values[k])
            # Add total to first axes
            
            axes[0].step(galaxy_time[::-1],flow_values['all'][::-1], where='mid',alpha=0.4,color=line_dict[args.model[i]])
            #all_median = medfilt(flow_values['all'][::-1],5)
            bin_medians, bin_edges, binnumber = stats.binned_statistic(galaxy_time[::-1],flow_values['all'][::-1],statistic='median',bins=bins)
            flow_arrays[args.model[i]] = np.ma.array(bin_medians, mask=np.isnan(bin_medians)),
            print(bin_medians)
            axes[0].plot(0.5*(bins[1:]+bins[:-1]),np.ma.array(bin_medians, mask=np.isnan(bin_medians)),marker='o', markersize=4, label=args.model[i],color=line_dict[args.model[i]])

            # Now compute phase contributions and add to specific plots
            hot = 100*(flow_values['hot'][::-1]/flow_values['all'][::-1])
            warm_ionised = 100*((flow_values['hot'][::-1]+flow_values['warm_ionised'][::-1])/flow_values['all'][::-1])
            warm_neutral = 100*((flow_values['hot'][::-1]+flow_values['warm_ionised'][::-1]+flow_values['warm_neutral'][::-1])/flow_values['all'][::-1])
            cold = 100*((flow_values['hot'][::-1]+flow_values['warm_ionised'][::-1]+flow_values['warm_neutral'][::-1]+flow_values['cold'][::-1])/flow_values['all'][::-1])
            
            print(args.model[i], cold[0],warm_neutral[0],warm_ionised[0],hot[0])
            axes[2+i].plot(galaxy_time[::-1],cold,marker='o', markersize=2, label=phase_labels['cold'],color=phase_colours['cold'])
            axes[2+i].fill_between(galaxy_time[::-1],cold,color=phase_colours['cold'],alpha=0.2)
            
            axes[2+i].plot(galaxy_time[::-1],warm_neutral,marker='o', markersize=2, label=phase_labels['warm_neutral'],color=phase_colours['warm_neutral'])
            axes[2+i].fill_between(galaxy_time[::-1],warm_neutral,color=phase_colours['warm_neutral'],alpha=0.2)
            
            axes[2+i].plot(galaxy_time[::-1],warm_ionised,marker='o', markersize=2, label=phase_labels['warm_ionised'],color=phase_colours['warm_ionised'])
            axes[2+i].fill_between(galaxy_time[::-1],warm_ionised,color=phase_colours['warm_ionised'],alpha=0.2)

            axes[2+i].plot(galaxy_time[::-1],hot,marker='o', markersize=2, label=phase_labels['hot'],color=phase_colours['hot'])
            axes[2+i].fill_between(galaxy_time[::-1],hot,color=phase_colours['hot'],alpha=0.2)
            axes[2+i].text(0.65, 0.3, args.model[i],
                        transform=axes[2+i].transAxes, fontsize=16,verticalalignment='top',
                        color='black')
            axes[2+i].tick_params(labelsize=14)
            axes[2+i].xaxis.set_ticks_position('both')
            axes[2+i].yaxis.set_ticks_position('both')
            axes[2+i].minorticks_on()
            axes[2+i].tick_params(which='both',axis="both",direction="in",labelbottom='off')

        # Add legends to axes
        ax = axes[0]
        #ax.set_ylim([0.01,30])
        ax.set_ylabel(r'$\dot{M}_{\rm %s}$ [$M_{\odot}$/yr]'%args.flowtype, fontsize=16)
        ax.tick_params(labelsize=14)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
        ax.set_yscale('log')
        ax.legend(loc='best',fontsize=16,frameon=False)
        ax.text(0.8, 0.95, r'$r= %.1f R_{\rm vir,DM}$'%(float(args.r)),
                        transform=ax.transAxes, fontsize=16,verticalalignment='top',
                        color='black')
        ax = axes[1]
        ax.set_ylabel(r'$\dot{M}_{\rm %s}/\dot{M}_{\rm %s,%s}$'%(args.flowtype,args.flowtype,args.model[-1]), fontsize=16)
        ax.tick_params(labelsize=14)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
        #ax.set_yscale('log')
        # Add ratio to second axes
        for i in range(0, len(args.model)-1):
            ratio = flow_arrays[args.model[i]][0]/flow_arrays[args.model[-1]][0]
            ax.plot(0.5*(bins[1:]+bins[:-1]),ratio,marker='o', markersize=4, label='Ratio with %s'%(args.model[i]),color=line_dict[args.model[i]])
        ax.legend(loc='best',fontsize=16,frameon=False)

        # Add top ticks for redshift
        ax = axes[0]
        axR = ax.twiny()
        maxt = cosmo.age(args.maxz).value
        mint = cosmo.age(13.0).value
        ax.set_xlim(mint, maxt)
        axR.set_xlim(mint, maxt)
        topticks1 = np.array([1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
        topticks1 = topticks1[topticks1 >= args.maxz]
        topticks2 = cosmo.age(topticks1).value
        axR.set_xticklabels(topticks1)
        axR.set_xticks(topticks2)
        axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
        axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
        axR.set_xlabel(r'$z$', fontsize=18)
        axR.tick_params(labelsize=14)

        ax = axes[2]
        ax.legend(loc='best',fontsize=16,frameon=False,ncol=2)

        # Add label for rest of axes
        ax = axes[-1]
        ax.set_xlabel(r'$t$ [Gyr]', fontsize=18)
        ax.tick_params(which='both',axis="both",direction="in",bottom=True)
        fig.text(0.01, 0.25, r'Mass fraction [$\%$]', fontsize=18,
                rotation='vertical',color='black',va='center')

        if args.NUT:
            args.ind = 'NUT'
        fig.subplots_adjust(top=0.95,bottom=0.05,left=0.1,right=0.98)
        fig.savefig(os.getcwd()+'/instantaneous_ratio_'+args.flowtype+'_'+args.var+'_'+str(args.ind)+'_'+d_key+'rvir.png', format='png', dpi=200)
    
    elif args.type == 'outflow_escape':
        from unyt import G
        nrow = int(1+len(args.model))
        height_ratios = []
        height_ratios.append(0.5)
        for i in range(0, len(args.model)):
            height_ratios.append(0.3)
        figsize = plt.figaspect(float((5.0 * len(args.model)) / (5.0 * 2)))
        fig = plt.figure(figsize=2*figsize, facecolor='w', edgecolor='k')
        plot_grid = fig.add_gridspec(nrow, 1, wspace=0, hspace=0, height_ratios=height_ratios)
        axes = []
        axes.append(fig.add_subplot(plot_grid[0]))
        for i in range(1,nrow):
            axes.append(fig.add_subplot(plot_grid[i],sharex=axes[0]))
        axes = np.asarray(axes)

        support_vars = {'grav_therpfrsphere':'r','grav_magpfrsphere':'m','grav_crpfrsphere':'forestgreen','grav_totpfrsphere':'black'}
        support_lstyle = {'grav_therpfrsphere':'-','grav_magpfrsphere':'-','grav_crpfrsphere':'-','grav_totpfrsphere':':'}
        # Get data
        plt_setting = plotting_dictionary['v_sphere_r']
        for i in range(0, len(args.model)):
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
                args.model[i] = args.model[i].replace('_','\_')
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
                args.model[i] = args.model[i].replace('_','\_')
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            files = os.listdir(groupspath)
            ozyfiles = []

            for f in files:
                if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                    ozyfiles.append(f)
            
            ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
            galaxy_time = []
            galaxy_masses = []
            flow_values = {'all':[],'hot':[],'warm_ionised':[],'warm_neutral':[],'cold':[]}
            support_values = {'all':{'grav_therpfrsphere':[],'grav_magpfrsphere':[],'grav_crpfrsphere':[],'grav_totpfrsphere':[]},
                                'hot':{'grav_therpfrsphere':[],'grav_magpfrsphere':[],'grav_crpfrsphere':[],'grav_totpfrsphere':[]},
                                'warm_ionised':{'grav_therpfrsphere':[],'grav_magpfrsphere':[],'grav_crpfrsphere':[],'grav_totpfrsphere':[]},
                                'warm_neutral':{'grav_therpfrsphere':[],'grav_magpfrsphere':[],'grav_crpfrsphere':[],'grav_totpfrsphere':[]},
                                'cold':{'grav_therpfrsphere':[],'grav_magpfrsphere':[],'grav_crpfrsphere':[],'grav_totpfrsphere':[]}}
            v_escape = []

            progind = args.ind
            if args.NUT:
                sim = ozy.load(os.path.join(groupspath, ozyfiles[0]))
                virial_mass = [i.virial_quantities['mass'] for i in sim.galaxies]
                progind = np.argmax(virial_mass)
            
            for ozyfile in ozyfiles:
                if progind == -1:
                    continue
                # Load OZY file
                try:
                    sim = ozy.load(os.path.join(groupspath, ozyfile))
                except:
                    continue
                # Initialise simulation parameters
                redshift = sim.simulation.redshift
                if redshift <= args.maxz:
                    try:
                        progind = sim.galaxies[progind].progen_galaxy_star
                    except:
                        progind = -1
                    continue
                h = sim.simulation.hubble_constant
                cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

                # Age of Universe at this redshift
                thubble = cosmo.age(redshift).value

                galaxy_time.append(thubble)

                gal = sim.galaxies[progind]
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
                    try:
                        flow_values[-2] = flow_values[-1]
                    except:
                        pass
                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    progind = -1
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.5*args.dr,'rvir'),
                                    rmax=(args.r+0.5*args.dr,'rvir'),save=True,recompute=args.recompute,
                                    remove_subs=args.rm_subs,pdf_bins=args.nbins)
                
                d_key = str(int(100*args.r))
                try:
                    for ftype in flow_values.keys():
                        d = gf.data['v_sphere_r_'+d_key+'rvir_'+ftype].in_units(plt_setting['units'])
                        if args.flowtype == 'inflow':
                            d = -d
                        flow_values[ftype].append(d[2])
                        for sv in support_vars.keys():
                            try:
                                d = gf.data[sv+'_'+d_key+'rvir_'+ftype]
                                support_values[ftype][sv].append(d[2])
                            except:
                                pass
                except:
                    raise KeyError('This variable is not present in this GalacticFlow object: '+str(gf.data[args.var+'_'+d_key+'rvir']))

                # Add escape velocity
                ves = 2*G*gal.mass['total']/(0.2*gal.halo.virial_quantities['radius'])
                ves = np.sqrt(ves).to('km/s').d
                print('V escape: ',ves, 0.2*gal.halo.virial_quantities['radius'].to('kpc'),args.model[i],redshift)
                v_escape.append(ves)
            
            # Convert to numpy arrays
            for k in flow_values.keys():
                flow_values[k] = np.asarray(flow_values[k])
                for s in support_vars.keys():
                    support_values[k][s] = np.asarray(support_values[k][s])

            # Add total to first axes
            axes[0].plot(galaxy_time[::-1],flow_values[args.var][::-1,1],marker='o', markersize=4, label=names[args.model[i]],color=line_dict[args.model[i]])
            axes[0].fill_between(galaxy_time[::-1],flow_values[args.var][::-1,3],flow_values[args.var][::-1,4],alpha=0.4,color=line_dict[args.model[i]])
            
            # Add escape velocity to first axes
            axes[0].plot(galaxy_time[::-1],v_escape[::-1],linestyle='--',color=line_dict[args.model[i]])

            # Now add support for each pressure
            for sv in support_vars.keys():
                plot_settings = plotting_dictionary[sv]
                try:
                    axes[1+i].plot(galaxy_time[::-1],support_values[args.var][sv][::-1,0],marker='o', markersize=2, 
                                    label=plot_settings['label'],color=support_vars[sv],linestyle=support_lstyle[sv])
                except:
                    pass
            axes[1+i].text(0.65, 0.3, names[args.model[i]],
                        transform=axes[1+i].transAxes, fontsize=16,verticalalignment='top',
                        color='black')
            axes[1+i].tick_params(labelsize=14)
            axes[1+i].xaxis.set_ticks_position('both')
            axes[1+i].yaxis.set_ticks_position('both')
            axes[1+i].minorticks_on()
            axes[1+i].tick_params(which='both',axis="both",direction="in",labelbottom='off')
            axes[1+i].set_yscale('symlog')
            axes[1+i].set_ylim([-1e3,4e2])

            # Add equilibrium line for pressure support
            maxt = cosmo.age(args.maxz).value
            mint = cosmo.age(13.0).value
            axes[1+i].fill_between([mint,maxt],[-1,-1],[1,1],alpha=0.4,color='b')

        # Add legends to axes
        ax = axes[0]
        # ax.set_ylim([2e-2,50])
        ax.set_ylabel(plt_setting['label'], fontsize=16)
        ax.tick_params(labelsize=14)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
        ax.set_yscale('log')
        ax.legend(loc='best',fontsize=16,frameon=False)
        ax.text(0.8, 0.95, r'$r= %.1f R_{\rm vir,DM}$'%(float(args.r)),
                        transform=ax.transAxes, fontsize=16,verticalalignment='top',
                        color='black')

        ax = axes[-1]
        ax.legend(loc='upper left',fontsize=16,frameon=False,ncol=4)

        # Add top ticks for redshift
        ax = axes[0]
        axR = ax.twiny()
        maxt = cosmo.age(args.maxz).value
        mint = cosmo.age(13.0).value
        ax.set_xlim(mint, maxt)
        axR.set_xlim(mint, maxt)
        topticks1 = np.array([1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
        topticks1 = topticks1[topticks1 >= args.maxz]
        topticks2 = cosmo.age(topticks1).value
        axR.set_xticklabels(topticks1)
        axR.set_xticks(topticks2)
        axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
        axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
        axR.set_xlabel(r'$z$', fontsize=18)
        axR.tick_params(labelsize=14)

        # Add label for rest of axes
        ax = axes[-1]
        ax.set_xlabel(r'$t$ [Gyr]', fontsize=18)
        ax.tick_params(which='both',axis="both",direction="in",bottom=True)
        fig.text(0.01, 0.25, r'Radial pressure support', fontsize=18,
                rotation='vertical',color='black',va='center')

        if args.NUT:
            args.ind = 'NUT'
        fig.subplots_adjust(top=0.95,bottom=0.05,left=0.1,right=0.98)
        if args.rm_subs:
            fig.savefig(os.getcwd()+'/outflow_escape_'+args.flowtype+'_'+args.var+'_'+str(args.ind)+'_'+d_key+'rvir_rmsubs.png', format='png', dpi=200)
        else:
            fig.savefig(os.getcwd()+'/outflow_escape_'+args.flowtype+'_'+args.var+'_'+str(args.ind)+'_'+d_key+'rvir.png', format='png', dpi=200)
    
    elif args.type == 'outflow_escape_phases':
        from unyt import G
        fig, axes = plt.subplots(2,2, figsize=(13,8),dpi=300,facecolor='w',edgecolor='k',sharex=True,sharey=True)

        # Get data
        plt_setting = plotting_dictionary['v_sphere_r']
        for i in range(0, len(args.model)):
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
                args.model[i] = args.model[i].replace('_','\_')
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
                args.model[i] = args.model[i].replace('_','\_')
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            files = os.listdir(groupspath)
            ozyfiles = []

            for f in files:
                if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                    ozyfiles.append(f)
            
            ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
            galaxy_time = []
            galaxy_masses = []
            flow_values = {'hot':[],'warm_ionised':[],'warm_neutral':[],'cold':[]}
            v_escape = []

            progind = args.ind
            if args.NUT:
                sim = ozy.load(os.path.join(groupspath, ozyfiles[0]))
                virial_mass = [i.virial_quantities['mass'] for i in sim.galaxies]
                progind = np.argmax(virial_mass)
            
            for ozyfile in ozyfiles:
                if progind == -1:
                    continue
                # Load OZY file
                try:
                    sim = ozy.load(os.path.join(groupspath, ozyfile))
                except:
                    continue
                # Initialise simulation parameters
                redshift = sim.simulation.redshift
                if redshift <= args.maxz:
                    try:
                        progind = sim.galaxies[progind].progen_galaxy_star
                    except:
                        progind = -1
                    continue
                h = sim.simulation.hubble_constant
                cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

                # Age of Universe at this redshift
                thubble = cosmo.age(redshift).value

                galaxy_time.append(thubble)

                gal = sim.galaxies[progind]
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
                    for ftype in flow_values.keys():
                        flow_values[ftype][-2] = flow_values[ftype][-1]
                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    progind = -1
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.5*args.dr,'rvir'),
                                    rmax=(args.r+0.5*args.dr,'rvir'),save=True,recompute=args.recompute,
                                    remove_subs=args.rm_subs,pdf_bins=args.nbins)
                
                d_key = str(int(100*args.r))
                try:
                    for ftype in flow_values.keys():
                        d = gf.data['v_sphere_r_'+d_key+'rvir_'+ftype].in_units(plt_setting['units'])
                        if args.flowtype == 'inflow':
                            d = -d
                        flow_values[ftype].append(d[2])
                except:
                    raise KeyError('This variable is not present in this GalacticFlow object: '+str(gf.data['v_sphere_r_'+d_key+'rvir']))

                # Add escape velocity
                ves = 2*G*gal.mass['total']/(0.2*gal.halo.virial_quantities['radius'])
                ves = np.sqrt(ves).to('km/s').d
                print('V escape: ',ves, 0.2*gal.halo.virial_quantities['radius'].to('kpc'),args.model[i],redshift)
                if bad_mass == True:
                    print('Changing ves')
                    v_escape[-1] = v_escape[-2]
                v_escape.append(ves)
            
            # Convert to numpy arrays
            for k in flow_values.keys():
                flow_values[k] = np.asarray(flow_values[k])

            for k,key in enumerate(flow_values.keys()):
                ax = axes[int(k/2),k%2]

                # Add gas velocity to each phase axes
                print(flow_values[key][::-1,1])
                ax.plot(galaxy_time[::-1],medfilt(flow_values[key][::-1,1],5),linewidth=2, label=names[args.model[i]],color=line_dict[args.model[i]])
                ax.fill_between(galaxy_time[::-1],medfilt(flow_values[key][::-1,3],5),medfilt(flow_values[key][::-1,4],5),alpha=0.2,color=line_dict[args.model[i]])
                
                # Add escape velocity to each phase axes
                ax.plot(galaxy_time[::-1],v_escape[::-1],linestyle='--',linewidth=2,color=line_dict[args.model[i]])


        # Add legend of simulation types
        ax = axes[0,0]
        first_legend = ax.legend(loc='lower right', fontsize=14,frameon=False,ncol=len(args.model))
        ax.add_artist(first_legend)
        # Add legend of velocity type (escape and gas outflow)
        ax = axes[0,1]
        dummy_lines = []

        dummy_lines.append(ax.plot([],[], color="black", ls = '-',label = r'$v_{\rm out}$')[0])
        dummy_lines.append(ax.plot([],[], color="black", ls = '--',label = r'$v_{\rm esc}$')[0])
        second_legend = ax.legend(handles=dummy_lines, loc='lower left', frameon=False, fontsize=14,ncol=2)         
        ax.add_artist(second_legend)

        p = 0
        for i in range(0,2):
            for j in range(0,2):
                ax = axes[i,j]
                ax.tick_params(labelsize=14)
                ax.xaxis.set_ticks_position('both')
                ax.yaxis.set_ticks_position('both')
                ax.minorticks_on()
                ax.tick_params(which='both',axis="both",direction="in")
                ax.set_yscale('log')
                phase = list(flow_values.keys())[p]
                ax.text(0.02, 0.95, r'\textbf{%s}'%phase_labels[phase],
                        transform=ax.transAxes, fontsize=20,verticalalignment='top',
                        color=phase_colours[phase], weight='bold')
                ax.set_ylim([10,450])
                p += 1

        axes[0,0].set_ylabel(plt_setting['label'], fontsize=20)
        axes[1,0].set_ylabel(plt_setting['label'], fontsize=20)
        
        axes[0,1].text(0.75, 0.95, r'$r= %.1f R_{\rm vir,DM}$'%(float(args.r)),
                        transform=axes[0,1].transAxes, fontsize=16,verticalalignment='top',
                        color='black')

        # Add top ticks for redshift
        for i in range(0,2):
            ax = axes[0,i]
            axR = ax.twiny()
            maxt = cosmo.age(args.maxz).value
            mint = cosmo.age(13.0).value
            ax.set_xlim(mint, maxt)
            axR.set_xlim(mint, maxt)
            topticks1 = np.array([1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
            topticks1 = topticks1[topticks1 >= args.maxz]
            topticks2 = cosmo.age(topticks1).value
            axR.set_xticklabels(topticks1)
            axR.set_xticks(topticks2)
            axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
            axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
            axR.set_xlabel(r'$z$', fontsize=20)
            axR.tick_params(labelsize=16)

        # Add label for rest of axes
        for i in range(0,2):
            ax = axes[1,i]
            ax.set_xlabel(r'$t$ [Gyr]', fontsize=20)
            ax.tick_params(which='both',axis="both",direction="in",bottom=True,labelsize=16)

        if args.NUT:
            args.ind = 'NUT'
        fig.subplots_adjust(top=0.92,bottom=0.1,left=0.06,right=0.98,hspace=0,wspace=0)
        if args.rm_subs:
            fig.savefig(os.getcwd()+'/outflow_escape_phases_'+args.flowtype+'_'+str(args.ind)+'_'+d_key+'rvir_rmsubs.pdf', format='pdf', dpi=300)
        else:
            fig.savefig(os.getcwd()+'/outflow_escape_phases_'+args.flowtype+'_'+str(args.ind)+'_'+d_key+'rvir.pdf', format='pdf', dpi=300)
    elif args.type == 'outflowrate_escape_phases':
        from unyt import G
        fig, axes = plt.subplots(4,2, figsize=(13,11),dpi=300,facecolor='w',edgecolor='k',sharex=True,sharey='row')

        # Get data
        plt_setting1 = plotting_dictionary['massflow_rate']
        plt_setting2 = plotting_dictionary['v_sphere_r']
        for i in range(0, len(args.model)):
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
                args.model[i] = args.model[i].replace('_','\_')
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
                args.model[i] = args.model[i].replace('_','\_')
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            files = os.listdir(groupspath)
            ozyfiles = []

            for f in files:
                if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                    ozyfiles.append(f)
            
            ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
            galaxy_time = []
            galaxy_masses = []
            flow_values = {'hot':[[],[]],'warm_ionised':[[],[]],'warm_neutral':[[],[]],'cold':[[],[]]}
            v_escape = []

            progind = args.ind
            if args.NUT:
                sim = ozy.load(os.path.join(groupspath, ozyfiles[0]))
                virial_mass = [i.virial_quantities['mass'] for i in sim.galaxies]
                progind = np.argmax(virial_mass)
            
            for ozyfile in ozyfiles:
                if progind == -1:
                    continue
                # Load OZY file
                try:
                    sim = ozy.load(os.path.join(groupspath, ozyfile))
                except:
                    continue
                # Initialise simulation parameters
                redshift = sim.simulation.redshift
                if redshift <= args.maxz:
                    try:
                        progind = sim.galaxies[progind].progen_galaxy_star
                    except:
                        progind = -1
                    continue
                h = sim.simulation.hubble_constant
                cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

                # Age of Universe at this redshift
                thubble = cosmo.age(redshift).value

                galaxy_time.append(thubble)

                gal = sim.galaxies[progind]
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
                    for ftype in flow_values.keys():
                        flow_values[ftype][0][-2] = flow_values[ftype][0][-1]
                        flow_values[ftype][1][-2] = flow_values[ftype][1][-1]
                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    progind = -1
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.5*args.dr,'rvir'),
                                    rmax=(args.r+0.5*args.dr,'rvir'),save=True,recompute=args.recompute,
                                    remove_subs=args.rm_subs,pdf_bins=args.nbins)
                
                d_key = str(int(100*args.r))
                try:
                    for ftype in flow_values.keys():
                        d = gf.data['massflow_rate_'+d_key+'rvir_'+ftype].in_units(plt_setting1['units'])
                        if args.flowtype == 'inflow':
                            d = -d
                        flow_values[ftype][0].append(d)
                        d = gf.data['v_sphere_r_'+d_key+'rvir_'+ftype].in_units(plt_setting2['units'])
                        if args.flowtype == 'inflow':
                            d = -d
                        flow_values[ftype][1].append(d[2])
                except:
                    raise KeyError('This variable is not present in this GalacticFlow object: '+str(gf.data['v_sphere_r_'+d_key+'rvir']))

                # Add escape velocity
                ves = 2*G*gal.mass['total']/(0.2*gal.halo.virial_quantities['radius'])
                ves = np.sqrt(ves).to('km/s').d
                if bad_mass == True:
                    print('Changing ves')
                    v_escape[-1] = v_escape[-2]
                v_escape.append(ves)
            
            # Convert to numpy arrays
            for k in flow_values.keys():
                flow_values[k][0] = np.asarray(flow_values[k][0])
                flow_values[k][1] = np.asarray(flow_values[k][1])

            for k,key in enumerate(flow_values.keys()):

                # 1. Add masflow rates
                ax = axes[2*int(k/2),k%2]

                # Add gas velocity to each phase axes
                ax.step(galaxy_time[::-1],flow_values[key][0][::-1], where='mid',alpha=0.4,color=line_dict[args.model[i]])
                all_median = medfilt(flow_values[key][0][::-1],5)
                ax.plot(galaxy_time[::-1],all_median, linewidth=2,label=names[args.model[i]],color=line_dict[args.model[i]])

                # 2. Add gas velocities
                ax = axes[2*int(k/2)+1,k%2]

                # Add gas velocity to each phase axes
                ax.plot(galaxy_time[::-1],medfilt(flow_values[key][1][::-1,1],5),linewidth=2, label=names[args.model[i]],color=line_dict[args.model[i]])
                ax.fill_between(galaxy_time[::-1],medfilt(flow_values[key][1][::-1,3],5),medfilt(flow_values[key][1][::-1,4],5),alpha=0.2,color=line_dict[args.model[i]])
                
                # Add escape velocity to each phase axes
                ax.plot(galaxy_time[::-1],v_escape[::-1],linestyle='--',linewidth=2,color=line_dict[args.model[i]])


        # Add legend of simulation types
        ax = axes[0,0]
        first_legend = ax.legend(loc='lower right', fontsize=14,frameon=False,ncol=len(args.model))
        ax.add_artist(first_legend)
        # Add legend of velocity type (escape and gas outflow)
        ax = axes[1,0]
        dummy_lines = []

        dummy_lines.append(ax.plot([],[], color="black", ls = '-',label = r'$v_{\rm outflow}$')[0])
        dummy_lines.append(ax.plot([],[], color="black", ls = '--',label = r'$v_{\rm esc}$')[0])
        second_legend = ax.legend(handles=dummy_lines, loc='lower left', frameon=False, fontsize=14,ncol=2)         
        ax.add_artist(second_legend)

        p = 0
        for i in range(0,2):
            for j in range(0,2):
                ax = axes[2*i,j]
                ax.tick_params(labelsize=14)
                ax.xaxis.set_ticks_position('both')
                ax.yaxis.set_ticks_position('both')
                ax.minorticks_on()
                ax.tick_params(which='both',axis="both",direction="in")
                ax.set_yscale('log')
                phase = list(flow_values.keys())[p]
                ax.text(0.02, 0.95, r'\textbf{%s}'%phase_labels[phase],
                        transform=ax.transAxes, fontsize=20,verticalalignment='top',
                        color=phase_colours[phase], weight='bold')
                ax.set_ylim([5e-3,20])
                ax = axes[2*i+1,j]
                ax.tick_params(labelsize=14)
                ax.xaxis.set_ticks_position('both')
                ax.yaxis.set_ticks_position('both')
                ax.minorticks_on()
                ax.tick_params(which='both',axis="both",direction="in")
                ax.set_yscale('log')
                phase = list(flow_values.keys())[p]
                ax.text(0.02, 0.95, r'\textbf{%s}'%phase_labels[phase],
                        transform=ax.transAxes, fontsize=20,verticalalignment='top',
                        color=phase_colours[phase], weight='bold')
                ax.set_ylim([25,450])
                p += 1
                

        axes[0,0].set_ylabel(r'$\dot{M}_{\rm outflow}$ [M$_{\odot}$/yr]', fontsize=20)
        axes[1,0].set_ylabel(plt_setting2['label'], fontsize=20)
        axes[2,0].set_ylabel(r'$\dot{M}_{\rm outflow}$ [M$_{\odot}$/yr]', fontsize=20)
        axes[3,0].set_ylabel(plt_setting2['label'], fontsize=20)
        
        axes[0,1].text(0.75, 0.95, r'$r= %.1f R_{\rm vir,DM}$'%(float(args.r)),
                        transform=axes[0,1].transAxes, fontsize=16,verticalalignment='top',
                        color='black')

        # Add top ticks for redshift
        for i in range(0,2):
            ax = axes[0,i]
            axR = ax.twiny()
            maxt = cosmo.age(args.maxz).value
            mint = cosmo.age(13.0).value
            ax.set_xlim(mint, maxt)
            axR.set_xlim(mint, maxt)
            topticks1 = np.array([1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
            topticks1 = topticks1[topticks1 >= args.maxz]
            topticks2 = cosmo.age(topticks1).value
            axR.set_xticklabels(topticks1)
            axR.set_xticks(topticks2)
            axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
            axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
            axR.set_xlabel(r'$z$', fontsize=20)
            axR.tick_params(labelsize=16)

        # Add label for rest of axes
        for i in range(0,2):
            ax = axes[3,i]
            ax.set_xlabel(r'$t$ [Gyr]', fontsize=20)
            ax.tick_params(which='both',axis="both",direction="in",bottom=True,labelsize=16)

        if args.NUT:
            args.ind = 'NUT'
        fig.subplots_adjust(top=0.95,bottom=0.07,left=0.06,right=0.99,hspace=0,wspace=0)
        if args.rm_subs:
            fig.savefig(os.getcwd()+'/outflowrate_escape_phases_'+args.flowtype+'_'+str(args.ind)+'_'+d_key+'rvir_rmsubs.pdf', format='pdf', dpi=300)
        else:
            fig.savefig(os.getcwd()+'/outflowrate_escape_phases_'+args.flowtype+'_'+str(args.ind)+'_'+d_key+'rvir.pdf', format='pdf', dpi=300)
    elif args.type == 'outflowfraction_escape_phases':
        from unyt import G
        fig, axes = plt.subplots(4,2, figsize=(13,11),dpi=300,facecolor='w',edgecolor='k',sharex=True,sharey='row')

        # Get data
        plt_setting1 = plotting_dictionary['massflow_rate']
        plt_setting2 = plotting_dictionary['v_sphere_r']
        for i in range(0, len(args.model)):
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
                args.model[i] = args.model[i].replace('_','\_')
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
                args.model[i] = args.model[i].replace('_','\_')
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            files = os.listdir(groupspath)
            ozyfiles = []

            for f in files:
                if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                    ozyfiles.append(f)
            
            ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
            galaxy_time = []
            galaxy_masses = []
            flow_values = {'hot':[[],[]],'warm_ionised':[[],[]],'warm_neutral':[[],[]],'cold':[[],[]]}
            v_escape = []

            progind = args.ind
            if args.NUT:
                sim = ozy.load(os.path.join(groupspath, ozyfiles[0]))
                virial_mass = [i.virial_quantities['mass'] for i in sim.galaxies]
                progind = np.argmax(virial_mass)
            
            for ozyfile in ozyfiles:
                if progind == -1:
                    continue
                # Load OZY file
                try:
                    sim = ozy.load(os.path.join(groupspath, ozyfile))
                except:
                    continue
                # Initialise simulation parameters
                redshift = sim.simulation.redshift
                if redshift <= args.maxz:
                    try:
                        progind = sim.galaxies[progind].progen_galaxy_star
                    except:
                        progind = -1
                    continue
                print(args.model[i],ozyfile)
                h = sim.simulation.hubble_constant
                cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

                # Age of Universe at this redshift
                thubble = cosmo.age(redshift).value

                galaxy_time.append(thubble)

                gal = sim.galaxies[progind]
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
                    for ftype in flow_values.keys():
                        flow_values[ftype][0][-2] = flow_values[ftype][0][-1]
                        flow_values[ftype][1][-2] = flow_values[ftype][1][-1]
                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    progind = -1
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.5*args.dr,'rvir'),
                                    rmax=(args.r+0.5*args.dr,'rvir'),save=True,recompute=args.recompute,
                                    remove_subs=args.rm_subs,pdf_bins=args.nbins)
                
                d_key = str(int(100*args.r))
                try:
                    for ftype in flow_values.keys():
                        d = gf.data['massflow_rate_'+d_key+'rvir_'+ftype].in_units(plt_setting1['units'])
                        d = d / gf.data['massflow_rate_'+d_key+'rvir_all'].in_units(plt_setting1['units'])
                        if args.flowtype == 'inflow':
                            d = -d
                        flow_values[ftype][0].append(d)
                        d = gf.data['v_sphere_r_'+d_key+'rvir_'+ftype].in_units(plt_setting2['units'])
                        if args.flowtype == 'inflow':
                            d = -d
                        flow_values[ftype][1].append(d[2])
                except:
                    raise KeyError('This variable is not present in this GalacticFlow object: '+str(gf.data['v_sphere_r_'+d_key+'rvir']))

                # Add escape velocity
                ves = 2*G*gal.mass['total']/(0.2*gal.halo.virial_quantities['radius'])
                ves = np.sqrt(ves).to('km/s').d
                if bad_mass == True:
                    print('Changing ves')
                    v_escape[-1] = v_escape[-2]
                v_escape.append(ves)
            
            # Convert to numpy arrays
            for k in flow_values.keys():
                flow_values[k][0] = np.asarray(flow_values[k][0])
                flow_values[k][1] = np.asarray(flow_values[k][1])

            for k,key in enumerate(flow_values.keys()):

                # 1. Add masflow rates
                ax = axes[2*int(k/2),k%2]

                # Add gas velocity to each phase axes
                ax.step(galaxy_time[::-1],flow_values[key][0][::-1], where='mid',alpha=0.4,color=line_dict[args.model[i]])
                all_median = medfilt(flow_values[key][0][::-1],5)
                ax.plot(galaxy_time[::-1],all_median, linewidth=2,label=names[args.model[i]],color=line_dict[args.model[i]])

                # 2. Add gas velocities
                ax = axes[2*int(k/2)+1,k%2]

                # Add gas velocity to each phase axes
                ax.plot(galaxy_time[::-1],medfilt(flow_values[key][1][::-1,1],5),linewidth=2, label=names[args.model[i]],color=line_dict[args.model[i]])
                ax.fill_between(galaxy_time[::-1],medfilt(flow_values[key][1][::-1,3],5),medfilt(flow_values[key][1][::-1,4],5),alpha=0.2,color=line_dict[args.model[i]])
                
                # Add escape velocity to each phase axes
                ax.plot(galaxy_time[::-1],v_escape[::-1],linestyle='--',linewidth=2,color=line_dict[args.model[i]])


        # Add legend of simulation types
        ax = axes[0,0]
        first_legend = ax.legend(loc='lower right', fontsize=14,frameon=False,ncol=len(args.model))
        ax.add_artist(first_legend)
        # Add legend of velocity type (escape and gas outflow)
        ax = axes[1,0]
        dummy_lines = []

        dummy_lines.append(ax.plot([],[], color="black", ls = '-',label = r'$v_{\rm outflow}$')[0])
        dummy_lines.append(ax.plot([],[], color="black", ls = '--',label = r'$v_{\rm esc}$')[0])
        second_legend = ax.legend(handles=dummy_lines, loc='lower left', frameon=False, fontsize=14,ncol=2)         
        ax.add_artist(second_legend)

        p = 0
        for i in range(0,2):
            for j in range(0,2):
                ax = axes[2*i,j]
                ax.tick_params(labelsize=14)
                ax.xaxis.set_ticks_position('both')
                ax.yaxis.set_ticks_position('both')
                ax.minorticks_on()
                ax.tick_params(which='both',axis="both",direction="in")
                ax.set_yscale('log')
                phase = list(flow_values.keys())[p]
                ax.text(0.02, 0.95, r'\textbf{%s}'%phase_labels[phase],
                        transform=ax.transAxes, fontsize=20,verticalalignment='top',
                        color=phase_colours[phase], weight='bold')
                ax.set_ylim([5e-3,20])
                ax = axes[2*i+1,j]
                ax.tick_params(labelsize=14)
                ax.xaxis.set_ticks_position('both')
                ax.yaxis.set_ticks_position('both')
                ax.minorticks_on()
                ax.tick_params(which='both',axis="both",direction="in")
                ax.set_yscale('log')
                phase = list(flow_values.keys())[p]
                ax.text(0.02, 0.95, r'\textbf{%s}'%phase_labels[phase],
                        transform=ax.transAxes, fontsize=20,verticalalignment='top',
                        color=phase_colours[phase], weight='bold')
                ax.set_ylim([25,450])
                p += 1
                

        axes[0,0].set_ylabel(r'$\dot{\chi}_{\rm outflow}$', fontsize=20)
        axes[1,0].set_ylabel(plt_setting2['label'], fontsize=20)
        axes[2,0].set_ylabel(r'$\dot{\chi}_{\rm outflow}$', fontsize=20)
        axes[3,0].set_ylabel(plt_setting2['label'], fontsize=20)
        
        axes[0,1].text(0.75, 0.95, r'$r= %.1f R_{\rm vir,DM}$'%(float(args.r)),
                        transform=axes[0,1].transAxes, fontsize=16,verticalalignment='top',
                        color='black')

        # Add top ticks for redshift
        for i in range(0,2):
            ax = axes[0,i]
            axR = ax.twiny()
            maxt = cosmo.age(args.maxz).value
            mint = cosmo.age(13.0).value
            ax.set_xlim(mint, maxt)
            axR.set_xlim(mint, maxt)
            topticks1 = np.array([1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
            topticks1 = topticks1[topticks1 >= args.maxz]
            topticks2 = cosmo.age(topticks1).value
            axR.set_xticklabels(topticks1)
            axR.set_xticks(topticks2)
            axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
            axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
            axR.set_xlabel(r'$z$', fontsize=20)
            axR.tick_params(labelsize=16)

        # Add label for rest of axes
        for i in range(0,2):
            ax = axes[3,i]
            ax.set_xlabel(r'$t$ [Gyr]', fontsize=20)
            ax.tick_params(which='both',axis="both",direction="in",bottom=True,labelsize=16)

        if args.NUT:
            args.ind = 'NUT'
        fig.subplots_adjust(top=0.95,bottom=0.07,left=0.06,right=0.99,hspace=0,wspace=0)
        if args.rm_subs:
            fig.savefig(os.getcwd()+'/outflowfraction_escape_phases_'+args.flowtype+'_'+str(args.ind)+'_'+d_key+'rvir_rmsubs.pdf', format='pdf', dpi=300)
        else:
            fig.savefig(os.getcwd()+'/outflowfraction_escape_phases_'+args.flowtype+'_'+str(args.ind)+'_'+d_key+'rvir.pdf', format='pdf', dpi=300)
    elif args.type == 'outflow_support_phases':
        from unyt import G
        fig, axes = plt.subplots(3,2, figsize=(13,10),dpi=300,facecolor='w',edgecolor='k',sharex=True,sharey=True)

        
        support_vars = {'grav_therpfrsphere':'r','grav_magpfrsphere':'m','grav_crpfrsphere':'forestgreen','grav_totpfrsphere':'black'}
        support_lstyle = {'grav_therpfrsphere':'--','grav_magpfrsphere':'-.','grav_crpfrsphere':'-','grav_totpfrsphere':':'}
        # Get data
        plt_setting = plotting_dictionary['v_sphere_r']
        for i in range(0, len(args.model)):
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
                args.model[i] = args.model[i].replace('_','\_')
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
                args.model[i] = args.model[i].replace('_','\_')
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            files = os.listdir(groupspath)
            ozyfiles = []

            for f in files:
                if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                    ozyfiles.append(f)
            
            ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
            galaxy_time = []
            galaxy_masses = []
            flow_values = {'hot':[],'warm_ionised':[],'warm_neutral':[],'cold':[]}
            support_values = {'hot':{'grav_therpfrsphere':[],'grav_crpfrsphere':[]},
                                'warm_ionised':{'grav_therpfrsphere':[],'grav_crpfrsphere':[]},
                                'warm_neutral':{'grav_therpfrsphere':[],'grav_crpfrsphere':[]},
                                'cold':{'grav_therpfrsphere':[],'grav_crpfrsphere':[]}}

            progind = args.ind
            if args.NUT:
                sim = ozy.load(os.path.join(groupspath, ozyfiles[0]))
                virial_mass = [i.virial_quantities['mass'] for i in sim.galaxies]
                progind = np.argmax(virial_mass)
            
            for ozyfile in ozyfiles:
                if progind == -1:
                    continue
                # Load OZY file
                try:
                    sim = ozy.load(os.path.join(groupspath, ozyfile))
                except:
                    continue
                # Initialise simulation parameters
                redshift = sim.simulation.redshift
                if redshift <= args.maxz:
                    try:
                        progind = sim.galaxies[progind].progen_galaxy_star
                    except:
                        progind = -1
                    continue
                h = sim.simulation.hubble_constant
                cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

                # Age of Universe at this redshift
                thubble = cosmo.age(redshift).value

                galaxy_time.append(thubble)

                gal = sim.galaxies[progind]
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
                    try:
                        flow_values[-2] = flow_values[-1]
                    except:
                        pass
                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    progind = -1
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.5*args.dr,'rvir'),
                                    rmax=(args.r+0.5*args.dr,'rvir'),save=True,recompute=args.recompute,
                                    remove_subs=args.rm_subs,pdf_bins=args.nbins)
                
                d_key = str(int(100*args.r))
                try:
                    for ftype in flow_values.keys():
                        for sv in support_vars.keys():
                            try:
                                d = gf.data[sv+'_'+d_key+'rvir_'+ftype]
                                support_values[ftype][sv].append(d[2])
                            except:
                                pass
                except:
                    raise KeyError('This variable is not present in this GalacticFlow object: '+str(gf.data['v_sphere_r_'+d_key+'rvir']))

            
            # Convert to numpy arrays
            for k in flow_values.keys():
                for s in support_vars.keys():
                    support_values[k][s] = np.asarray(support_values[k][s])

            for k,key in enumerate(flow_values.keys()):

                if k == 0:
                    ax = axes[i,0]
                else:
                    ax = axes[k-1,1]
                    
                # Now add support for each pressure
                for sv in support_vars.keys():
                    plot_settings = plotting_dictionary[sv]
                    try:
                        ax.plot(galaxy_time[::-1],medfilt(support_values[key][sv][::-1,0],5),marker='o', markersize=4, 
                                        color=line_dict[args.model[i]],linestyle=support_lstyle[sv],
                                        markerfacecolor=support_vars[sv])
                        ax.fill_between(galaxy_time[::-1],medfilt(support_values[key][sv][::-1,3],5),
                                        medfilt(support_values[key][sv][::-1,4],5),alpha=0.2,
                                        color=line_dict[args.model[i]])
                
                    except:
                        pass


        # Add legend of simulation types
        ax = axes[0,0]
        dummy_lines = []
        for m in args.model:
            dummy_lines.append(ax.plot([],[], color=line_dict[m], ls = '-',marker='o',label=names[m])[0])
        first_legend = ax.legend(handles=dummy_lines, loc='lower left', frameon=False, fontsize=14,ncol=len(args.model)) 
        ax.add_artist(first_legend)

        # Add legend of support type
        ax = axes[0,1]
        dummy_lines = []
        for sv in support_vars.keys():
            plot_settings = plotting_dictionary[sv]
            dummy_lines.append(ax.plot([],[], color='black', ls = support_lstyle[sv],marker='o',label=plot_settings['label'],markerfacecolor=support_vars[sv])[0])
        second_legend = ax.legend(handles=dummy_lines, loc='lower right', frameon=False, fontsize=14,ncol=2)     
        ax.add_artist(second_legend)

        for i in range(0,3):
            for j in range(0,2):
                ax = axes[i,j]
                ax.tick_params(labelsize=14)
                ax.xaxis.set_ticks_position('both')
                ax.yaxis.set_ticks_position('both')
                ax.minorticks_on()
                ax.tick_params(which='both',axis="both",direction="in")
                ax.set_yscale('symlog')
                if j == 1:
                    p = i + 1
                    phase = list(flow_values.keys())[p]
                    ax.text(0.02, 0.95, phase_labels[phase],
                        transform=ax.transAxes, fontsize=16,verticalalignment='top',
                        color=phase_colours[phase])
                # Add equilibrium line for pressure support
                maxt = cosmo.age(args.maxz).value
                mint = cosmo.age(13.0).value
                ax.fill_between([mint,maxt],[-1,-1],[1,1],alpha=0.4,color='grey')

        axes[0,0].text(0.02, 0.95, phase_labels['hot'],
                        transform=axes[0,0].transAxes, fontsize=16,verticalalignment='top',
                        color=phase_colours['hot'])
        axes[1,0].set_ylabel(r'Radial pressure support', fontsize=18)
        axes[1,0].text(0.02, 0.95, phase_labels['hot'],
                        transform=axes[1,0].transAxes, fontsize=16,verticalalignment='top',
                        color=phase_colours['hot'])
        axes[2,0].text(0.02, 0.95, phase_labels['hot'],
                        transform=axes[2,0].transAxes, fontsize=16,verticalalignment='top',
                        color=phase_colours['hot'])

        ax.text(0.02, 0.95, phase_labels[phase],
                        transform=ax.transAxes, fontsize=16,verticalalignment='top',
                        color=phase_colours[phase])
        
        axes[0,1].text(0.75, 0.95, r'$r= %.1f R_{\rm vir,DM}$'%(float(args.r)),
                        transform=axes[0,1].transAxes, fontsize=16,verticalalignment='top',
                        color='black')

        # Add top ticks for redshift
        for i in range(0,2):
            ax = axes[0,i]
            axR = ax.twiny()
            maxt = cosmo.age(args.maxz).value
            mint = cosmo.age(13.0).value
            ax.set_xlim(mint, maxt)
            axR.set_xlim(mint, maxt)
            topticks1 = np.array([1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
            topticks1 = topticks1[topticks1 >= args.maxz]
            topticks2 = cosmo.age(topticks1).value
            axR.set_xticklabels(topticks1)
            axR.set_xticks(topticks2)
            axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
            axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
            axR.set_xlabel(r'$z$', fontsize=18)
            axR.tick_params(labelsize=14)

        # Add label for rest of axes
        for i in range(0,2):
            ax = axes[1,i]
            ax.set_xlabel(r'$t$ [Gyr]', fontsize=18)
            ax.tick_params(which='both',axis="both",direction="in",bottom=True)

        if args.NUT:
            args.ind = 'NUT'
        fig.subplots_adjust(top=0.92,bottom=0.07,left=0.07,right=0.98,hspace=0,wspace=0)
        if args.rm_subs:
            fig.savefig(os.getcwd()+'/outflow_support_phases_'+args.flowtype+'_'+str(args.ind)+'_'+d_key+'rvir_rmsubs.png', format='png', dpi=300)
        else:
            fig.savefig(os.getcwd()+'/outflow_support_phases_'+args.flowtype+'_'+str(args.ind)+'_'+d_key+'rvir.png', format='png', dpi=300)
    elif args.type == 'outflow_support_pos_phases':
        from unyt import G
        fig, axes = plt.subplots(3,2, figsize=(13,10),dpi=300,facecolor='w',edgecolor='k',sharex=True,sharey=True)

        support_vars = {'grav_therpfrspherepos':'r','grav_crpfrspherepos':'forestgreen'}
        support_lstyle = {'grav_therpfrspherepos':'--','grav_crpfrspherepos':'-'}
        # Get data
        for i in range(0, len(args.model)):
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
                args.model[i] = args.model[i].replace('_','\_')
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
                args.model[i] = args.model[i].replace('_','\_')
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            files = os.listdir(groupspath)
            ozyfiles = []

            for f in files:
                if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                    ozyfiles.append(f)
            
            ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
            galaxy_time = []
            galaxy_masses = []
            flow_values = {'hot':[],'warm_ionised':[],'warm_neutral':[],'cold':[]}
            support_values = {'hot':{'grav_therpfrspherepos':[],'grav_crpfrspherepos':[]},
                                'warm_ionised':{'grav_therpfrspherepos':[],'grav_crpfrspherepos':[]},
                                'warm_neutral':{'grav_therpfrspherepos':[],'grav_crpfrspherepos':[]},
                                'cold':{'grav_therpfrspherepos':[],'grav_crpfrspherepos':[]}}

            progind = args.ind
            if args.NUT:
                sim = ozy.load(os.path.join(groupspath, ozyfiles[0]))
                virial_mass = [i.virial_quantities['mass'] for i in sim.galaxies]
                progind = np.argmax(virial_mass)
            
            for ozyfile in ozyfiles:
                if progind == -1:
                    continue
                # Load OZY file
                try:
                    sim = ozy.load(os.path.join(groupspath, ozyfile))
                except:
                    continue
                # Initialise simulation parameters
                redshift = sim.simulation.redshift
                if redshift <= args.maxz:
                    try:
                        progind = sim.galaxies[progind].progen_galaxy_star
                    except:
                        progind = -1
                    continue
                h = sim.simulation.hubble_constant
                cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

                # Age of Universe at this redshift
                thubble = cosmo.age(redshift).value

                galaxy_time.append(thubble)

                gal = sim.galaxies[progind]
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
                    try:
                        for ftype in flow_values.keys():
                            for sv in support_vars.keys():
                                support_values[ftype][sv][-1] = support_values[ftype][sv][-2]
                    except:
                        pass
                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    progind = -1
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.5*args.dr,'rvir'),
                                    rmax=(args.r+0.5*args.dr,'rvir'),save=True,recompute=args.recompute,
                                    remove_subs=args.rm_subs,pdf_bins=args.nbins)
                
                d_key = str(int(100*args.r))
                try:
                    for ftype in flow_values.keys():
                        for sv in support_vars.keys():
                            try:
                                d = gf.data[sv+'_'+d_key+'rvir_'+ftype]
                                if d[2][0] == d[2][1] and d[2][2] == 0:
                                    support_values[ftype][sv].append(np.zeros(7))
                                else:
                                    support_values[ftype][sv].append(d[2])
                            except:
                                print('Failed in: ',ftype,sv,args.model[i],ozyfile)
                                pass
                except:
                    raise KeyError('This variable is not present in this GalacticFlow object: '+str(gf.data['v_sphere_r_'+d_key+'rvir']))

            
            # Convert to numpy arrays
            for k in flow_values.keys():
                for s in support_vars.keys():
                    support_values[k][s] = np.asarray(support_values[k][s])

            for k,key in enumerate(flow_values.keys()):

                if k == 0:
                    ax = axes[i,0]
                else:
                    ax = axes[k-1,1]
                    
                # Now add support for each pressure
                for sv in support_vars.keys():
                    plot_settings = plotting_dictionary[sv]
                    try:
                        ax.plot(galaxy_time[::-1],medfilt(support_values[key][sv][::-1,1],5),marker='o', markersize=5, 
                                        color=line_dict[args.model[i]],linestyle=support_lstyle[sv],
                                        markerfacecolor=support_vars[sv],linewidth=2)
                        ax.fill_between(galaxy_time[::-1],medfilt(support_values[key][sv][::-1,3],5),
                                        medfilt(support_values[key][sv][::-1,4],5),alpha=0.2,
                                        color=line_dict[args.model[i]])
                    except:
                        pass


        # Add legend of simulation types
        ax = axes[0,0]
        dummy_lines = []
        for m in args.model:
            dummy_lines.append(ax.plot([],[], color=line_dict[m], ls = '-',marker='o',label=names[m])[0])
        first_legend = ax.legend(handles=dummy_lines, loc='upper right', frameon=False, fontsize=14,ncol=len(args.model)) 
        ax.add_artist(first_legend)

        # Add legend of support type
        ax = axes[0,1]
        dummy_lines = []
        for sv in support_vars.keys():
            plot_settings = plotting_dictionary[sv]
            dummy_lines.append(ax.plot([],[], color='black', ls = support_lstyle[sv],marker='o',label=plot_settings['label'],markerfacecolor=support_vars[sv])[0])
        second_legend = ax.legend(handles=dummy_lines, loc='best', frameon=False, fontsize=16,ncol=2)     
        ax.add_artist(second_legend)

        for i in range(0,3):
            for j in range(0,2):
                ax = axes[i,j]
                ax.tick_params(labelsize=14)
                ax.xaxis.set_ticks_position('both')
                ax.yaxis.set_ticks_position('both')
                ax.minorticks_on()
                ax.tick_params(which='both',axis="both",direction="in")
                ax.set_yscale('symlog',linthresh=0.5,linscale=1)
                ax.set_ylim([0,15])
                ax.tick_params(labelsize=16)
                if j == 1:
                    p = i + 1
                    phase = list(flow_values.keys())[p]
                    ax.text(0.02, 0.95, r'\textbf{%s}'%phase_labels[phase],
                        transform=ax.transAxes, fontsize=20,verticalalignment='top',
                        color=phase_colours[phase],weight='bold')
                # Add equilibrium line for pressure support
                maxt = cosmo.age(args.maxz).value
                mint = cosmo.age(13.0).value
                ax.fill_between([mint,maxt],[0,0],[1,1],alpha=0.4,color='grey')
                
        axes[0,0].set_ylabel(r'$-\nabla_{r}P_x/\rho g_{r}$', fontsize=20)
        axes[0,0].text(0.02, 0.95, r'\textbf{%s}'%phase_labels['hot'],
                        transform=axes[0,0].transAxes, fontsize=20,verticalalignment='top',
                        color=phase_colours['hot'],weight='bold')
        axes[1,0].set_ylabel(r'$-\nabla_{r}P_x/\rho g_{r}$', fontsize=20)
        axes[1,0].text(0.02, 0.95, r'\textbf{%s}'%phase_labels['hot'],
                        transform=axes[1,0].transAxes, fontsize=20,verticalalignment='top',
                        color=phase_colours['hot'],weight='bold')
        axes[2,0].set_ylabel(r'$-\nabla_{r}P_x/\rho g_{r}$', fontsize=20)
        axes[2,0].text(0.02, 0.95, r'\textbf{%s}'%phase_labels['hot'],
                        transform=axes[2,0].transAxes, fontsize=20,verticalalignment='top',
                        color=phase_colours['hot'],weight='bold')
        
        axes[0,1].text(0.70, 0.68, r'$r= %.1f R_{\rm vir,DM}$'%(float(args.r)),
                        transform=axes[0,1].transAxes, fontsize=20,verticalalignment='top',
                        color='black')

        # Add top ticks for redshift
        for i in range(0,2):
            ax = axes[0,i]
            axR = ax.twiny()
            maxt = cosmo.age(args.maxz).value
            mint = cosmo.age(13.0).value
            ax.set_xlim(mint, maxt)
            axR.set_xlim(mint, maxt)
            topticks1 = np.array([1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
            topticks1 = topticks1[topticks1 >= args.maxz]
            topticks2 = cosmo.age(topticks1).value
            axR.set_xticklabels(topticks1)
            axR.set_xticks(topticks2)
            axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
            axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
            axR.set_xlabel(r'$z$', fontsize=20)
            axR.tick_params(labelsize=16)

        # Add label for rest of axes
        for i in range(0,2):
            ax = axes[2,i]
            ax.set_xlabel(r'$t$ [Gyr]', fontsize=18)
            ax.minorticks_on()
            ax.tick_params(which='both',axis="both",direction="in",bottom=True,labelsize=16)

        if args.NUT:
            args.ind = 'NUT'
        fig.subplots_adjust(top=0.92,bottom=0.07,left=0.07,right=0.98,hspace=0,wspace=0)
        if args.rm_subs:
            fig.savefig(os.getcwd()+'/outflow_support_pos_phases_'+args.flowtype+'_'+str(args.ind)+'_'+d_key+'rvir_rmsubs.pdf', format='pdf', dpi=300)
        else:
            fig.savefig(os.getcwd()+'/outflow_support_pos_phases_'+args.flowtype+'_'+str(args.ind)+'_'+d_key+'rvir.pdf', format='pdf', dpi=300)
    elif args.type == 'instantaneous_vr':
        nrow = 3
        fig = plt.figure(figsize=(9,8), facecolor='w', edgecolor='k')
        plot_grid = fig.add_gridspec(nrow, 1, wspace=0, hspace=0)
        axes = []
        axes.append(fig.add_subplot(plot_grid[0]))
        for i in range(1,nrow):
            axes.append(fig.add_subplot(plot_grid[i],sharex=axes[0]))
        axes = np.asarray(axes)

        # Get data
        for i in range(0, len(args.model)):
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
                args.model[i] = args.model[i].replace('_','\_')
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
                args.model[i] = args.model[i].replace('_','\_')
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            files = os.listdir(groupspath)
            ozyfiles = []

            for f in files:
                if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                    ozyfiles.append(f)
            
            ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
            galaxy_time = []
            galaxy_masses = []
            flow_values = []
            vr = []
            metallicity = []

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
                print(args.model[i],redshift,ozyfile)
                if redshift <= args.maxz:
                    try:
                        progind = sim.galaxies[progind].progen_galaxy_star
                    except:
                        progind = -1
                    continue
                h = sim.simulation.hubble_constant
                cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

                # Age of Universe at this redshift
                thubble = cosmo.age(redshift).value

                galaxy_time.append(thubble)

                gal = sim.galaxies[progind]
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
                    try:
                        flow_values[-2] = flow_values[-1]
                    except:
                        pass
                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    progind = -1
                
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.5*args.dr,'rvir'),
                                    rmax=(args.r+0.5*args.dr,'rvir'),save=True,recompute=args.recompute,
                                    remove_subs=args.rm_subs,pdf_bins=args.nbins)
                
                d_key = str(int(100*args.r))
                mrate = gf.data['massflow_rate_'+d_key+'rvir_all'].in_units('Msun/yr')
                vspeed = gf.data['v_sphere_r_'+d_key+'rvir_all'].in_units('km/s')
                metal = gf.data['metallicity_'+d_key+'rvir_all']
                if args.flowtype == 'inflow':
                    mrate = -mrate
                    vspeed = -vspeed
                flow_values.append(mrate)
                # We only take the outflow rate weighted value
                vr.append(vspeed[2])
                print(metal)
                metallicity.append(metal[2])
            
            # Convert to numpy arrays
            flow_values = np.asarray(flow_values)
            vr = np.asarray(vr)
            metallicity = np.asarray(metallicity)
            
            # Add total to first axes
            axes[0].step(galaxy_time[::-1],flow_values[::-1], where='mid',alpha=0.4,color=line_dict[args.model[i]])
            all_median = medfilt(flow_values[::-1],5)
            axes[0].plot(galaxy_time[::-1],all_median,marker='o', markersize=4, label=names[args.model[i]],color=line_dict[args.model[i]])
            
            # Add gas velocity
            axes[1].plot(galaxy_time[::-1],vr[::-1,1],marker='o', markersize=4, label=names[args.model[i]],color=line_dict[args.model[i]])
            axes[1].fill_between(galaxy_time[::-1],vr[::-1,3],vr[::-1,4],alpha=0.4,color=line_dict[args.model[i]])
            
            # Add gas metallicity
            axes[2].plot(galaxy_time[::-1],metallicity[::-1,1],marker='o', markersize=4, label=names[args.model[i]],color=line_dict[args.model[i]])
            axes[2].fill_between(galaxy_time[::-1],metallicity[::-1,3],metallicity[::-1,4],alpha=0.4,color=line_dict[args.model[i]])


        # Add legends to axes
        ax = axes[0]
        ax.set_ylim([1,50])
        ax.set_ylabel(r'$\dot{M}_{\rm %s}$ [$M_{\odot}$/yr]'%args.flowtype, fontsize=20)
        ax.tick_params(labelsize=14)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
        ax.set_yscale('log')
        ax.legend(loc='best',fontsize=16,frameon=False, ncol=len(args.model))
        ax.text(0.77, 0.95, r'$r= %.1f R_{\rm vir,DM}$'%(float(args.r)),
                        transform=ax.transAxes, fontsize=16,verticalalignment='top',
                        color='black')
        ax = axes[1]
        ax.set_ylabel(r'$v_{r}$ [km/s]', fontsize=20)
        ax.tick_params(labelsize=14)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
        ax.set_yscale('log')
        
        ax = axes[2]
        ax.set_ylabel(r'$Z$ [$Z_{\odot}$]', fontsize=16)
        ax.tick_params(labelsize=14)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
        ax.set_yscale('log')

        # Add top ticks for redshift
        ax = axes[0]
        axR = ax.twiny()
        maxt = cosmo.age(args.maxz).value
        mint = cosmo.age(13.0).value
        ax.set_xlim(mint, maxt)
        axR.set_xlim(mint, maxt)
        topticks1 = np.array([1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
        topticks1 = topticks1[topticks1 >= args.maxz]
        topticks2 = cosmo.age(topticks1).value
        axR.set_xticklabels(topticks1)
        axR.set_xticks(topticks2)
        axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
        axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
        axR.set_xlabel(r'$z$', fontsize=18)
        axR.tick_params(labelsize=14)

        # Add label for rest of axes
        ax = axes[-1]
        ax.set_xlabel(r'$t$ [Gyr]', fontsize=18)
        ax.tick_params(which='both',axis="both",direction="in",bottom=True)

        if args.NUT:
            args.ind = 'NUT'
        fig.subplots_adjust(top=0.93,bottom=0.09,left=0.1,right=0.98)
        if args.rm_subs:
            fig.savefig(os.getcwd()+'/instantaneous_vr_'+args.flowtype+'_'+str(args.ind)+'_'+d_key+'rvir_rmsubs.png', format='png', dpi=200)
        else:
            fig.savefig(os.getcwd()+'/instantaneous_vr_'+args.flowtype+'_'+str(args.ind)+'_'+d_key+'rvir.png', format='png', dpi=200)
    else:
        # Now add all to a single plot
        fig, ax = plt.subplots(1,1, figsize=(7,5),dpi=100,facecolor='w',edgecolor='k')
        plt_setting = plotting_dictionary[args.var]
        ax.set_xlabel(r'$t$ [Gyr]', fontsize=16)
        ax.set_ylabel(plt_setting['label'], fontsize=16)
        ax.tick_params(labelsize=12)
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.minorticks_on()
        ax.tick_params(which='both',axis="both",direction="in")
        ax.set_yscale('log')
        for i in range(0, len(args.model)):
            if args.model[i][0] != '/':
                simfolder = os.path.join(os.getcwd(), args.model[i])
            else:
                simfolder = args.model[i]
                args.model[i] = args.model[i].split('/')[-1]
            if not os.path.exists(simfolder):
                raise Exception('The given simulation name is not found in this directory!')
            
            groupspath = os.path.join(simfolder, 'Groups')
            files = os.listdir(groupspath)
            ozyfiles = []

            for f in files:
                if f.startswith('ozy_') and args.start < int(f[4:-5]) < args.end:
                    ozyfiles.append(f)
            
            ozyfiles = sorted(ozyfiles, key=lambda x:x[4:-5], reverse=True)
            galaxy_time = []
            galaxy_masses = []
            flow_values = []

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
                if redshift <= args.maxz:
                    try:
                        progind = sim.galaxies[progind].progen_galaxy_star
                    except:
                        progind = -1
                    continue
                h = sim.simulation.hubble_constant
                cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

                # Age of Universe at this redshift
                thubble = cosmo.age(redshift).value

                galaxy_time.append(thubble)

                gal = sim.galaxies[progind]
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
                    try:
                        flow_values[-2] = flow_values[-1]
                    except:
                        pass
                try:
                    progind = sim.galaxies[progind].progen_galaxy_star
                except:
                    progind = -1
                
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.5*args.dr,'rvir'),
                                    rmax=(args.r+0.5*args.dr,'rvir'),save=True,recompute=args.recompute,
                                    remove_subs=args.rm_subs,pdf_bins=args.nbins)
                
                d_key = str(int(100*args.r))
                try:
                    d = gf.data[args.var+'_'+d_key+'rvir_all'].in_units(plt_setting['units'])
                    if args.flowtype == 'inflow':
                        d = -d
                    flow_values.append(d)
                except:
                    raise KeyError('This variable is not present in this GalacticFlow object: '+str(gf.data[args.var+'_'+d_key+'rvir']))

            ax.plot(galaxy_time,flow_values,marker='o', markersize=2, label=args.model[i])

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
        ax.legend(loc='best', fontsize=14,frameon=False)

        if args.NUT:
            args.ind = 'NUT'
        fig.savefig(os.getcwd()+'/'+args.flowtype+'_'+args.var+'_'+str(args.ind)+'_'+d_key+'rvir.png', format='png', dpi=200)
            
