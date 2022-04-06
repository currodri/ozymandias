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
    parser.add_argument('--type', type=str, default='outflow', help='Flow variable to be plotted.')
    parser.add_argument('--ind', type=int, default=0, help='Index of the galaxy in the ending snapshot.')
    parser.add_argument('--var', type=str, default='massflow_rate', help='Flow variable to be plotted.')
    parser.add_argument('--r', type=float, default=0.2, help='Fraction of rvir at which flow measured.')
    parser.add_argument('--start', type=int, default=10, help='Starting index to look for galaxy.')
    parser.add_argument('--end', type=int, default=1000, help='Ending index to look for galaxy.')
    parser.add_argument('--maxz',type=float,default=1.5, help="Maximum redshift displayed in plot.")
    parser.add_argument('--flowtype', type=str, default='outflow', help='Flow type to be plotted.')
    parser.add_argument('--sfr', type=str, default='100Myr',help='What SFR indicator bin to use.')
    parser.add_argument('--NUT', type=bool, default=True, help='If True, it looks for NUT as the most massive galaxy (stars) in the last snapshot.')
    parser.add_argument('--recompute',action='store_true', help='If present, it recomputes 2D profiles.')

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

        line_dict = {'cosmoNUThd':'b','cosmoNUTmhd':'m','cosmoNUTcrmhd':'g','cosmoNUTcrmhd\_nost':'olive','cosmoNUTcrmhd\_noheat':'darkkhaki'}
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
                                                rmax=(0.5+0.01,'rvir'),save=False,recompute=args.recompute)
                gf['inflow_100'] = compute_flows(gal,os.path.join(groupspath, ozyfile),'inflow',rmin=(1.0-0.01,'rvir'),
                                                rmax=(1.0+0.01,'rvir'),save=False,recompute=args.recompute)
                gf['outflow_50'] = compute_flows(gal,os.path.join(groupspath, ozyfile),'outflow',rmin=(0.5-0.01,'rvir'),
                                                rmax=(0.5+0.01,'rvir'),save=False,recompute=args.recompute)
                gf['outflow_100'] = compute_flows(gal,os.path.join(groupspath, ozyfile),'outflow',rmin=(1.0-0.01,'rvir'),
                                                rmax=(1.0+0.01,'rvir'),save=False,recompute=args.recompute)
                
                
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
                
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.01,'rvir'),
                                    rmax=(args.r+0.01,'rvir'),save=True,recompute=args.recompute)
                
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

        line_dict = {'cosmoNUThd':'b','cosmoNUTmhd':'m','cosmoNUTcrmhd':'g','cosmoNUTcrmhd\_nost':'olive','cosmoNUTcrmhd\_noheat':'darkkhaki'}
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
                
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.01,'rvir'),
                                    rmax=(args.r+0.01,'rvir'),save=True,recompute=args.recompute,separate_phases=True)
                
                d_key = str(int(100*args.r))
                factor = 1
                for d in range(0,3):
                        for f in line_styles.keys():
                            v = plotting_vars[d][0]
                            plt_setting = plotting_dictionary[v]
                            print(v+'_'+d_key+'rvir_'+str(f))
                            try:
                                value = gf.data[v+'_'+d_key+'rvir_'+str(f)].in_units(plt_setting['units'])
                                if args.type == 'inflow':
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

        phase_labels = {'hot':'Hot','warm_ionised':'Warm ionised','warm_neutral':'Warm neutral','cold':'Cold'}
        phase_colours = {'hot':'r','warm_ionised':'orange','warm_neutral':'yellow','cold':'b'}
        line_dict = {'cosmoNUThd':'royalblue','cosmoNUTmhd':'m','cosmoNUTcrmhd':'g','cosmoNUTcrmhd\_nost':'olive','cosmoNUTcrmhd\_noheat':'darkgoldenrod'}
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
                
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.01,'rvir'),
                                    rmax=(args.r+0.01,'rvir'),save=True,recompute=args.recompute)
                
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

        phase_labels = {'hot':'Hot','warm_ionised':'Warm ionised','warm_neutral':'Warm neutral','cold':'Cold'}
        phase_colours = {'hot':'r','warm_ionised':'orange','warm_neutral':'yellow','cold':'b'}
        line_dict = {'cosmoNUThd':'royalblue','cosmoNUTmhd':'m','cosmoNUTcrmhd':'g','cosmoNUTcrmhd\_nost':'olive','cosmoNUTcrmhd\_noheat':'darkgoldenrod'}
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
                
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.01,'rvir'),
                                    rmax=(args.r+0.01,'rvir'),save=True,recompute=args.recompute)
                
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
                sfr.append(gal.sfr[args.sfr].in_units('Msun/yr'))
            
            # Convert to numpy arrays
            for k in flow_values.keys():
                flow_values[k] = np.asarray(flow_values[k])

            # Add total to first axes
            axes[0].step(galaxy_time[::-1],flow_values['all'][::-1], where='mid',alpha=0.4,color=line_dict[args.model[i]])
            all_median = medfilt(flow_values['all'][::-1],5)
            axes[0].plot(galaxy_time[::-1],all_median,marker='o', markersize=4, label=args.model[i],color=line_dict[args.model[i]])

            # Add SFR to second axes
            axes[1].plot(galaxy_time[::-1],sfr[::-1],marker='o', markersize=4, label=args.model[i],color=line_dict[args.model[i]])

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
        ax.set_ylabel(r'SFR [$M_{\odot}$/yr]', fontsize=16)
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
        topticks1 = np.array([1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 13.0])
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
        fig.savefig(os.getcwd()+'/instantaneous_sfr_'+args.flowtype+'_'+args.var+'_'+str(args.ind)+'_'+d_key+'rvir.png', format='png', dpi=200)
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

        phase_labels = {'hot':'Hot','warm_ionised':'Warm ionised','warm_neutral':'Warm neutral','cold':'Cold'}
        phase_colours = {'hot':'r','warm_ionised':'orange','warm_neutral':'yellow','cold':'b'}
        line_dict = {'cosmoNUThd':'royalblue','cosmoNUTmhd':'m','cosmoNUTcrmhd':'g','cosmoNUTcrmhd\_nost':'olive','cosmoNUTcrmhd\_noheat':'darkgoldenrod'}
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
                
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.flowtype,rmin=(args.r-0.01,'rvir'),
                                    rmax=(args.r+0.01,'rvir'),save=True,recompute=args.recompute)
                
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
        topticks1 = np.array([1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 13.0])
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
                
                gf = compute_flows(gal,os.path.join(groupspath, ozyfile),args.type,rmin=(args.r-0.01,'rvir'),
                                    rmax=(args.r+0.01,'rvir'),save=True,recompute=args.recompute)
                
                d_key = str(int(100*args.r))
                try:
                    d = gf.data[args.var+'_'+d_key+'rvir_all'].in_units(plt_setting['units'])
                    if args.type == 'inflow':
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
        fig.savefig(os.getcwd()+'/'+args.type+'_'+args.var+'_'+str(args.ind)+'_'+d_key+'rvir.png', format='png', dpi=200)
            