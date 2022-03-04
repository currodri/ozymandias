"""
This code allows the comparison of the halo and galaxy metallicity between different
simulations and redshifts.

By: F. Rodriguez Montero (15/02/2022)
"""

# Import required libraries
from numpy.core.fromnumeric import var
import ozy
import numpy as np
import os
import sys
import subprocess
import argparse
from astropy.cosmology import FlatLambdaCDM, z_at_value
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/amr')
from amr2 import io_ramses,filtering,amr_integrator
from ozy.utils import init_region
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
sns.set(style="white")
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# hfont = {'fontname':'Helvetica'}
# matplotlib.rc('text', usetex = True)
# matplotlib.rc('font', **{'family' : "serif"})
# params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
# matplotlib.rcParams.update(params)


if __name__ == '__main__':

    # Parse the command line arguments.
    parser = argparse.ArgumentParser(description='Tracking of galaxy metallicity across cosmic time.')
    parser.add_argument('model', type=str, nargs='+', help='Simulation names from which extract OZY results.')
    parser.add_argument('--ind', type=int, default=0, help='Index of the galaxy in the ending snapshot.')
    parser.add_argument('--zinit',type=float,default=2.0, help="Redshift at which to start the tracking.")
    parser.add_argument('--zend',type=float,default=6.0, help="Redshift at which to end the tracking.")
    parser.add_argument('--NUT', type=bool, default=True, help='If True, it looks for NUT as the most massive galaxy (stars) in the last snapshot.')
    args = parser.parse_args()

    if not isinstance(args.model, list):
        args.model = [args.model]


    # Now add all to a simple plot
    fig, axes = plt.subplots(1, 1, sharex=True, figsize=(6,5), dpi=100, facecolor='w', edgecolor='k')
    ax = axes
    ax.set_xlabel(r'$Z_{\rm gas,galaxy} [Z_{\odot}]$', fontsize=16)
    ax.set_ylabel(r'$Z_{\rm gas,halo} [Z_{\odot}]$', fontsize=16)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.tick_params(labelsize=12)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.minorticks_on()
    ax.tick_params(which='both',axis="both",direction="in")

    for i in range(0, len(args.model)):
        
        if args.model[i][0] != '/':
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
        for j in range(0, len(args.z)):
            os.chdir(simfolder)
            result = subprocess.run('IDtoZetas.out -ask '+str(int(args.z[j])), shell=True,stdout=subprocess.PIPE)
            os.chdir('../')
            indexout = int(result.stdout.decode('utf-8').split('is')[1].split('(z')[0])

            ozyfile = 'ozy_%05d.hdf5' % (indexout)
            sim = ozy.load(os.path.join(groupspath, ozyfile))

            progind = args.ind
            if args.NUT:
                virial_mass = [i.virial_quantities['mass'] for i in sim.galaxies]
                progind = np.argmax(virial_mass)

            gal = sim.galaxies[progind]

            # Get galaxy metallicity from ozy
            zgal = gal.metallicity['gas']

            # Compute halo metallicity
            output_path = sim.simulation.fullpath
            selected_reg = init_region(gal,'sphere', rmin=(0.2,'rvir'),rmax=(1.0,'rvir'))
            all_filt = filtering.filter()
            glob_attrs = amr_integrator.amr_region_attrs()
            glob_attrs.nvars = 1
            glob_attrs.nwvars = 1
            glob_attrs.nfilter = 1
            amr_integrator.allocate_amr_regions_attrs(glob_attrs)
            glob_attrs.varnames.T.view('S128')[0] = 'metallicity'.ljust(128)
            glob_attrs.wvarnames.T.view('S128')[0] = 'mass'.ljust(128)
            glob_attrs.filters[0] = all_filt

            # Begin integration
            amr_integrator.integrate_region(output_path,selected_reg,glob_attrs)
            zhalo = glob_attrs.data[0,0,0,0]

            if not sim.simulation.physics['cr'] and not sim.simulation.physics['magnetic']:
                ax.plot(zgal,zhalo, marker='s', markeredgecolor='b', markerfacecolor='none',label=args.model[i])
            elif not sim.simulation.physics['cr'] and sim.simulation.physics['magnetic']:
                ax.plot(zgal,zhalo, marker='s', markeredgecolor='m', markerfacecolor='none',label=args.model[i])
            elif sim.simulation.physics['cr']:
                ax.plot(zgal,zhalo, marker='s', markeredgecolor='g', markerfacecolor='none',label=args.model[i])
    fig.subplots_adjust(top=0.91, bottom=0.1,right=0.97,hspace=0.0)
    if args.NUT:
        args.ind = 'NUT'
    fig.savefig(os.getcwd()+'/comp_metals_'+str(args.ind)+'.png', format='png', dpi=300)