"""
This is a simple code that compares three different phase diagrams using a single plot.
These can be taken from:
    (a) the same redshift but different galaxies
    (b) the same galaxy but different redshifts
    (c) the same galaxy, same redshift but different simulations

By: F. Rodriguez Montero (23/02/2021)
"""

# Import required libraries
import ozy
import numpy as np
import os
import argparse
import subprocess
from ozy.loader import Galaxy
from unyt import G
import ozy
from ozy.phase_diagrams import compute_phase_diagram, \
                                plot_compare_phase_diagram, \
                                plot_compare_stacked_pd
from ozy.utils import closest_snap_z, find_neigh_snaps, get_tdyn
from astropy.cosmology import FlatLambdaCDM

# Name conversions
# TODO: This shouldn't be here...
names = {'cosmoNUThd':'HD',
        'cosmoNUThd\_all\_cr10':'HDcr10',
        'cosmoNUThd\_all\_cr20':'HDcr20',
        'cosmoNUTmhd':'MHD',
        'cosmoNUTcrmhd':'CRMHD',
        'cosmoNUTcrmhd\_nost':'nsCRMHD',
        'cosmoNUTcrmhd\_noheat':'nhCRMHD',
        'cosmoNUTrticrmhd':'RTCRiMHD',
        'cosmoNUTcrmhd\_3e29':'3e29CRMHD'
        }

if __name__ == '__main__':

    # Parse the command line arguments.
    parser = argparse.ArgumentParser(description='COMPARE PHASE DIAGRAMS')
    parser.add_argument('type', type=str, default='model', help='Comparison option.')
    parser.add_argument('--ind', type=int, default=0, nargs='+', help='Index of the galaxy at the desired redshift.')
    parser.add_argument('--z', type=float, default=[3.0], nargs='+', help='Redshift at which compare.')
    parser.add_argument('--model', type=str, nargs='+', help='Model names to compare.')
    parser.add_argument('--region',type=str, default='galaxy',help='Region from where to extract data.')
    parser.add_argument('--xvar', type=str, default='density', help='Horizontal axis variable.')
    parser.add_argument('--yvar', type=str, default='temperature', help='Vertical axis variable.')
    parser.add_argument('--field', type=str, default='gas/mass', help='Field used in the colormap.')
    parser.add_argument('--weight', type=str, default='gas/cumulative', help='Weighting variable.')
    parser.add_argument('--doflows',action='store_true', help='If present, it separates for the outflows and inflows.')
    parser.add_argument('--do_sf',action='store_true', help='If present, it separates the satr forming gas.')
    parser.add_argument('--sfeff', type=float, default=0.01, help='For the case of plots with SF eff, is the theshold included.')
    parser.add_argument('--layout', type=str, default='compact', help='How the plots should be organised.')
    parser.add_argument('--recompute',action='store_true', help='If present, it recomputes 2D profiles.')
    parser.add_argument('--scaletype',type=str, default='log_even', help='Type of scale for the bins.')
    parser.add_argument('--stats',type=str, default='none', help='What stats to use for the overplotted line.')
    parser.add_argument('--NUT', type=bool, default=True, help='If True, it looks for NUT as the most massive galaxy (stars) in the last snapshot.')
    args = parser.parse_args()


    if args.type == 'model':
        if not isinstance(args.model, list):
            print('You need multiple models if you want to compare models!')
            exit
        
        if args.doflows: print('You have asked for outflow/inflow contours!')
        if args.do_sf: print('You have asked for the contours of SF gas.')

        for z in range(0, len(args.z)):
            pds = []
            simfolders = []
            origfolder = os.getcwd()
            for i in range(0, len(args.model)):
            
                if args.model[i][0] != '/':
                    simfolder = os.path.join(os.getcwd(), args.model[i])
                    args.model[i] = args.model[i].replace('_','\_')
                else:
                    simfolder = args.model[i]
                    args.model[i] = args.model[i].split('/')[-1]
                    args.model[i] = args.model[i].replace('_','\_')
                simfolders.append(simfolder)
                #TODO: This should be done in a much cleaner way...
                print(args.model[i])
                if args.model[i] == 'cosmoNUTcrmhd': 
                    cr_flags = [True,True]
                elif args.model[i] == 'cosmoNUTcrmhd\_nost':
                    cr_flags = [False,False]
                elif args.model[i] == 'cosmoNUTcrmhd\_noheat':
                    cr_flags = [True,False]
                else:
                    cr_flags = [False,False]

                print(cr_flags)

                if not os.path.exists(simfolder):
                    raise Exception('The given simulation name is not found in this directory!: %s'%simfolder)
                
                groupspath = os.path.join(simfolder, 'Groups')
                os.chdir(simfolder)
                result = subprocess.run('IDtoZetas.out -ask '+str(args.z[z]), shell=True,stdout=subprocess.PIPE)
                os.chdir('../')
                indexout = int(result.stdout.decode('utf-8').split('is')[1].split('(z')[0])

                ozyfile = 'ozy_%05d.hdf5' % (indexout)

                sim = ozy.load(os.path.join(groupspath, ozyfile))

                progind = args.ind
                if args.NUT:
                    virial_mass = [i.virial_quantities['mass'] for i in sim.galaxies]
                    progind = np.argmax(virial_mass)
                    args.ind = 'NUT'
                gal = sim.galaxies[progind]

                if args.region == 'galaxy':
                    rmin, rmax = (0.0, 'rvir'),(0.2, 'rvir')
                elif args.region == 'halo':
                    rmin, rmax = (0.2, 'rvir'),(1.0, 'rvir')
                elif args.region == 'shell_rvir':
                    rmin, rmax = (0.99, 'rvir'),(1.1, 'rvir')
                elif args.region == 'shell_galaxy':
                    rmin, rmax = (0.19, 'rvir'),(0.21, 'rvir')
                print('scale type: ',args.scaletype)
                if args.doflows:
                    # Escape velocity at 0.2 Rvir of the halo
                    v_escape = 2*G*gal.mass['total']/(0.2*sim.halos[gal.parent_halo_index].virial_quantities['radius'])
                    v_escape = np.sqrt(v_escape).to('km/s').d
                    print('v_sphere_r/>=/%.3f/km*s**-1'%v_escape)
                    outflow = compute_phase_diagram(gal,os.path.join(groupspath, ozyfile), args.xvar,args.yvar, [args.field],
                                [args.weight],save=True,recompute=args.recompute, filter_conds='v_sphere_r/>/10/km*s**-1', filter_name='outflow'+'_'+args.region,
                                rmin=rmin,rmax=rmax,scaletype=args.scaletype,
                                cr_st = cr_flags[0],cr_heat=cr_flags[1])
                    inflow = compute_phase_diagram(gal,os.path.join(groupspath, ozyfile), args.xvar,args.yvar, [args.field],
                                [args.weight],save=True,recompute=args.recompute, filter_conds='v_sphere_r/</-10/km*s**-1', filter_name='inflow'+'_'+args.region,
                                rmin=rmin,rmax=rmax,scaletype=args.scaletype,
                                cr_st = cr_flags[0],cr_heat=cr_flags[1])
                    escape = compute_phase_diagram(gal,os.path.join(groupspath, ozyfile), args.xvar,args.yvar, [args.field],
                                [args.weight],save=True,recompute=args.recompute, filter_conds='v_sphere_r/>=/%.3f/km*s**-1'%v_escape, filter_name='escaping'+'_'+args.region,
                                rmin=rmin,rmax=rmax,scaletype=args.scaletype,
                                cr_st = cr_flags[0],cr_heat=cr_flags[1])
                if args.do_sf:
                    if sim.simulation.physics['magnetic'] or sim.simulation.physics['cr']:
                        sf_model = 'eff_FKmag'
                    else:
                        sf_model = 'eff_FK2'
                    condition_sf = sf_model+'/>=/%.3f/dimensionless'%args.sfeff
                    print(condition_sf)
                    starforming = compute_phase_diagram(gal,os.path.join(groupspath, ozyfile), args.xvar,args.yvar, [args.field],
                                    [args.weight],save=True,recompute=args.recompute, filter_conds=condition_sf, filter_name='starforming'+'_'+args.region,
                                    rmin=rmin,rmax=rmax,scaletype=args.scaletype,
                                    cr_st = cr_flags[0],cr_heat=cr_flags[1])
                pd = compute_phase_diagram(gal,os.path.join(groupspath, ozyfile), args.xvar,args.yvar, [args.field],
                            [args.weight],save=True,recompute=args.recompute,rmin=rmin,rmax=rmax,filter_name=args.region,scaletype=args.scaletype,
                            cr_st = cr_flags[0],cr_heat=cr_flags[1])
                if args.doflows and args.do_sf:
                    pds.append([pd,outflow,inflow,escape,starforming])
                elif args.doflows:
                    pds.append([pd,outflow,inflow,escape])
                elif args.do_sf:
                    pds.append([pd,starforming])
                else:
                    pds.append(pd)
                os.chdir(origfolder)
            namephase = '_'+args.xvar+'_'+args.yvar+'_'
            labels = [names[m] for m in args.model]
            plot_compare_phase_diagram(pds,args.field.split('/')[1],'compare_pd'+namephase+args.region+'_'+args.field.split('/')[1]+'_'+str(args.z[z]),
                                        weightvar=args.weight.split('/')[1],stats=args.stats,extra_labels=labels,
                                        doflows=args.doflows,do_sf=args.do_sf,layout=args.layout,gent=True)
            args.model = simfolders
    elif args.type == 'model_stacked':
        if not isinstance(args.model, list):
            print('You need multiple models if you want to compare models!')
            exit
        
        if args.doflows: print('You have asked for outflow/inflow contours!')
        if args.do_sf: print('You have asked for the contours of SF gas.')
        print('Phase diagrams are the result of stacking!')

        for z in range(0, len(args.z)):
            pds = []
            fields = []
            weights = []
            simfolders = []
            origfolder = os.getcwd()
            print('Running from %s'%origfolder)
            sample_tdyn = np.zeros(len(args.model))
            for i in range(0, len(args.model)):
                if args.model[i][0] != '/':
                    simfolder = os.path.join(os.getcwd(), args.model[i])
                    args.model[i] = args.model[i].replace('_','\_')
                else:
                    simfolder = args.model[i]
                    args.model[i] = args.model[i].split('/')[-1]
                    args.model[i] = args.model[i].replace('_','\_')
                simfolders.append(simfolder)

                if not os.path.exists(simfolder):
                    raise Exception('The given simulation name is not found in this directory!: %s'%simfolder)
                
                groupspath = os.path.join(simfolder, 'Groups')
                ozyfile = closest_snap_z(simfolder,args.z[z])

                sim = ozy.load(os.path.join(groupspath, ozyfile))

                progind = args.ind
                if args.NUT:
                    virial_mass = [i.virial_quantities['mass'] for i in sim.galaxies]
                    progind = np.argmax(virial_mass)
                    args.ind = 'NUT'
                gal = sim.galaxies[progind]

                tdyn = get_tdyn(gal).to('Gyr').d
                sample_tdyn[i] = tdyn
            max_tdyn = np.max(sample_tdyn)
            std_tdyn = np.std(sample_tdyn)
            print(sample_tdyn)
            redshift = sim.simulation.redshift
            h = sim.simulation.hubble_constant
            cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                                    Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)
            print('For z=%.1f the maximum dynamical time is %.3f Gyr with a STD of %.4f'%(args.z[z],max_tdyn,std_tdyn))
            print('The age of the Universe at z=%.1f is %.3f Gyr'%(args.z[z],cosmo.age(redshift).value))
            for i in range(0, len(args.model)):
                pds.append([])
                weights.append([])
                simfolder = simfolders[i]
                if not os.path.exists(simfolder):
                    raise Exception('The given simulation name is not found in this directory!: %s'%simfolder)
                
                groupspath = os.path.join(simfolder, 'Groups')
                ozyfile = closest_snap_z(simfolder,args.z[z])

                sim = ozy.load(os.path.join(groupspath, ozyfile))

                progind = args.ind
                if args.NUT:
                    virial_mass = [i.virial_quantities['mass'] for i in sim.galaxies]
                    progind = np.argmax(virial_mass)
                    args.ind = 'NUT'
                gal = sim.galaxies[progind]

                #TODO: This should be done in a much cleaner way...
                print(args.model[i])
                if args.model[i] == 'cosmoNUTcrmhd': 
                    cr_flags = [True,True]
                elif args.model[i] == 'cosmoNUTcrmhd\_nost':
                    cr_flags = [False,False]
                elif args.model[i] == 'cosmoNUTcrmhd\_noheat':
                    cr_flags = [True,False]
                else:
                    cr_flags = [False,False]

                print(cr_flags)

                orig_path = sim.simulation.fullpath.split('/')[-1]
                neigh_snaps, neigh_weights = find_neigh_snaps(simfolder, orig_path,
                                                                max_tdyn, returnweight=True)
                if args.region == 'galaxy':
                    rmin, rmax = (0.0, 'rvir'),(0.2, 'rvir')
                elif args.region == 'halo':
                    rmin, rmax = (0.2, 'rvir'),(1.0, 'rvir')
                elif args.region == 'shell_rvir':
                    rmin, rmax = (0.99, 'rvir'),(1.1, 'rvir')
                fields.append(args.field.split('/')[1])
                for n in range(0, len(neigh_snaps)):
                    sim = ozy.load(os.path.join(groupspath, neigh_snaps[n]))
                    if args.NUT:
                        virial_mass = [i.virial_quantities['mass'] for i in sim.galaxies]
                        progind = np.argmax(virial_mass)
                        args.ind = 'NUT'
                    gal = sim.galaxies[progind]
                    if args.doflows:
                        # Escape velocity at 0.2 Rvir of the halo
                        v_escape = 2*G*gal.mass['total']/(0.2*sim.halos[gal.parent_halo_index].virial_quantities['radius'])
                        v_escape = np.sqrt(v_escape).to('km/s').d
                        print('v_sphere_r/>=/%.3f/km*s**-1'%v_escape)
                        outflow = compute_phase_diagram(gal,os.path.join(groupspath, neigh_snaps[n]), args.xvar,args.yvar, [args.field],
                                    [args.weight],save=True,recompute=args.recompute, filter_conds='v_sphere_r/>/10/km*s**-1', filter_name='outflow'+'_'+args.region,
                                    rmin=rmin,rmax=rmax,scaletype=args.scaletype,
                                cr_st = cr_flags[0],cr_heat=cr_flags[1])
                        inflow = compute_phase_diagram(gal,os.path.join(groupspath, neigh_snaps[n]), args.xvar,args.yvar, [args.field],
                                    [args.weight],save=True,recompute=args.recompute, filter_conds='v_sphere_r/</-10/km*s**-1', filter_name='inflow'+'_'+args.region,
                                    rmin=rmin,rmax=rmax,scaletype=args.scaletype,
                                cr_st = cr_flags[0],cr_heat=cr_flags[1])
                        escape = compute_phase_diagram(gal,os.path.join(groupspath, neigh_snaps[n]), args.xvar,args.yvar, [args.field],
                                    [args.weight],save=True,recompute=args.recompute, filter_conds='v_sphere_r/>=/%.3f/km*s**-1'%v_escape, filter_name='escaping'+'_'+args.region,
                                    rmin=rmin,rmax=rmax,scaletype=args.scaletype,
                                cr_st = cr_flags[0],cr_heat=cr_flags[1])
                    if args.do_sf:
                        if sim.simulation.physics['magnetic'] or sim.simulation.physics['cr']:
                            sf_model = 'eff_FKmag'
                        else:
                            sf_model = 'eff_FK2'
                        condition_sf = sf_model+'/>=/%.3f/dimensionless'%args.sfeff
                        print(condition_sf)
                        starforming = compute_phase_diagram(gal,os.path.join(groupspath, neigh_snaps[n]), args.xvar,args.yvar, [args.field],
                                        [args.weight],save=True,recompute=args.recompute, filter_conds=condition_sf, filter_name='starforming'+'_'+args.region,
                                        rmin=rmin,rmax=rmax,scaletype=args.scaletype,
                                        cr_st = cr_flags[0],cr_heat=cr_flags[1])
                        
                    pd = compute_phase_diagram(gal,os.path.join(groupspath, neigh_snaps[n]), args.xvar,args.yvar, [args.field],
                                [args.weight],save=True,recompute=args.recompute,rmin=rmin,rmax=rmax,filter_name=args.region,scaletype=args.scaletype,
                                cr_st = cr_flags[0],cr_heat=cr_flags[1])
                    if args.doflows and args.do_sf:
                        pds[i].append([pd,outflow,inflow,escape,starforming])
                        weights[i].append([neigh_weights[n]]*5)
                    elif args.doflows:
                        pds[i].append([pd,outflow,inflow,escape])
                        weights[i].append([neigh_weights[n]]*4)
                    elif args.do_sf:
                        pds[i].append([pd,starforming])
                        weights[i].append([neigh_weights[n]]*2)
                    else:
                        pds[i].append(pd)
                        weights[i].append(neigh_weights[n])
                os.chdir(origfolder)

            print('Now plotting...')
            namephase = '_'+args.xvar+'_'+args.yvar+'_'
            labels = [names[m] for m in args.model]
            plot_compare_stacked_pd(pds,weights,fields,'compare_stacked_pd'+namephase+args.region+'_'+args.field.split('/')[1]+'_'+str(args.z[z]),
                                        weightvar=args.weight.split('/')[1],stats='mean',extra_labels=labels,
                                        doflows=args.doflows,do_sf=args.do_sf,layout=args.layout,gent=True)
            args.model = simfolders
    else:
        print('That comparison mode is not suported. Please check!')
        exit
        