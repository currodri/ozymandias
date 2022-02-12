"""
This is a simple code that allows the comparison of AMR and particle projections
between different snapshots/simulations.

By: F. Rodriguez Montero (24/02/2022)
"""

import os
import argparse
import subprocess
import numpy as np
import ozy
from ozy.utils import as_si
from ozy.projections import plot_comp_fe
if __name__ == '__main__':

    # Parse the command line arguments.
    parser = argparse.ArgumentParser(description='COMPARE PROJECTIONS')
    parser.add_argument('type', type=str, default='model_grid', help='Comparison option.')
    parser.add_argument('--ind', type=int, default=0, nargs='+', help='Index of the galaxy at the desired redshift.')
    parser.add_argument('--z', type=float, default=[2.0], nargs='+', help='Redshift at which compare.')
    parser.add_argument('--model', type=str, nargs='+', help='Model names to compare.')
    parser.add_argument('--fields', type=str, nargs='+', help='Fields used.')
    parser.add_argument('--NUT', type=bool, default=True, help='If True, it looks for NUT as the most massive galaxy (stars) in the last snapshot.')
    args = parser.parse_args()

    if args.type == 'model_grid':
        if not isinstance(args.model, list):
            print('You need multiple models if you want to compare models!')
            exit
        
        faceon_fits = []
        edgeon_fits = []
        labels = []

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

            if args.NUT:
                args.ind = 'NUT'
            faceon_fits.append(groupspath+'/%s_%05d_faceon.fits' % (args.ind,indexout))
            edgeon_fits.append(groupspath+'/%s_%05d_edgeon.fits' % (args.ind,indexout))
            # Get galaxy and model details for labels
            ozyfile = ozy.load(groupspath+'/ozy_%05d.hdf5'%indexout)
            if args.NUT:
                virial_mass = [i.virial_quantities['mass'] for i in ozyfile.galaxies]
                gal = ozyfile.galaxies[np.argmax(virial_mass)]
            stellar = r'$M_*$'+r'$={0:s}$'.format(as_si(gal.mass['stellar'].to('Msun').d,2)) + r' $M_{\odot}$'
            gas = r'$M_{\rm gas}$'+r'$={0:s}$'.format(as_si(gal.mass['gas'].to('Msun').d,2)) + r' $M_{\odot}$'
            labels.append([args.model[i],stellar,gas])
        fig = plot_comp_fe(faceon_fits,edgeon_fits,args.fields,returnfig=True,labels=labels)
        print('%s_comp_facevsedge_z%.2f.png' % (args.ind,args.z[0]))
        fig.savefig('%s_comp_facevsedge_z%.2f.png' % (args.ind,args.z[0]),dpi=300)