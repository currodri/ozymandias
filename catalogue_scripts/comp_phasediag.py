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
from yt import YTQuantity
import ozy
from ozy.phase_diagrams import compute_phase_diagram, plot_compare_phase_diagram,plot_single_phase_diagram

if __name__ == '__main__':

    # Parse the command line arguments.
    parser = argparse.ArgumentParser(description='COMPARE PHASE DIAGRAMS')
    parser.add_argument('type', type=str, default='model', help='Comparison option.')
    parser.add_argument('--ind', type=int, default=0, nargs='+', help='Index of the galaxy at the desired redshift.')
    parser.add_argument('--z', type=float, default=3.0, nargs='+', help='Redshift at which compare.')
    parser.add_argument('--model', type=str, nargs='+', help='Model names to compare.')
    parser.add_argument('--field', type=str, default='gas/mass', help='Field used in the colormap.')
    parser.add_argument('--weight', type=str, default='gas/cumulative', help='Weighting variable.')
    parser.add_argument('--NUT', type=bool, default=True, help='If True, it looks for NUT as the most massive galaxy (stars) in the last snapshot.')
    args = parser.parse_args()


    if args.type == 'model':
        if not isinstance(args.model, list):
            print('You need multiple models if you want to compare models!')
            exit
        
        pds = []

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
            os.chdir(simfolder)
            result = subprocess.run('IDtoZetas.out -ask '+str(args.z[0]), shell=True,stdout=subprocess.PIPE)
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
            pd = compute_phase_diagram(gal,os.path.join(groupspath, ozyfile), 'density','temperature', [args.field],
                        [args.weight],save=True,recompute=False)
            pds.append(pd)
        
        plot_compare_phase_diagram(pds,args.field.split('/')[1],'compare_pd_'+args.field.split('/')[1]+'_'+str(args.z[0])+'.png',weightvar=args.weight.split('/')[1],stats='mean',extra_labels=args.model)
    else:
        print('That comparison mode is not suported. Please check!')
        exit
        