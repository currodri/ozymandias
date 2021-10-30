"""
This code allows the comparison of the measured Stellar-to-Halo Mass Function in a cosmological simulation
to the prediction of models.

By: F. Rodriguez Montero (28/10/2021)
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
from yt import YTArray, YTQuantity
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/amr')
from amr2 import io_ramses
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


def Moster2013(z, mhalo):
    """
    This function returns the best fit results by Moster et al. (2013) for a given redshift z.

    Masses are returned in units of Msun.
    """

    # Taken from Table 1 of Moster et al. (2013)
    M10 = 11.590
    M10_err = 0.236
    M11 = 1.195
    M11_err = 0.353
    N10 = 0.0351
    N10_err = 0.0058
    N11 = -0.0247
    N11_err = 0.0069
    B10 = 1.376
    B10_err = 0.153
    B11 = -0.826
    B11_err = 0.225
    G10 = 0.608
    G10_err = 0.059
    G11 = 0.329
    G11_err = 0.173

    rf = z / (1+z)
    # Equation 11 in Moster et al. (2013)
    M1 = 10**(M10 + M11*rf)
    
    # Equation 12 in Moster et al. (2013)
    N = N10 + N11 * rf

    # Equation 13 in Moster et al. (2013)
    beta = B10 + B11 * rf

    # Equation 14 in Moster et al. (2013)

    gamma = G10 + G11 * rf

    # Equation 2 in Moster et al. (2013)
    m = 2 * N * ((mhalo / M1)**(-beta) + (mhalo / M1)**gamma)**(-1.0)
    m = m * mhalo

    return m


if __name__ == '__main__':

    # Parse the command line arguments.
    parser = argparse.ArgumentParser(description='Tracking of galaxy masses across cosmic time.')
    parser.add_argument('model', type=str, nargs='+', help='Simulation names from which extract OZY results.')
    parser.add_argument('--ind', type=int, default=0, help='Index of the galaxy in the ending snapshot.')
    parser.add_argument('--z',type=float,default=3.0, help="Redshift at which to extract data.")
    parser.add_argument('--ref', type=str, default='Moster+2013', help='Model to which compare.')
    parser.add_argument('--NUT', type=bool, default=True, help='If True, it looks for NUT as the most massive galaxy (stars) in the last snapshot.')
    args = parser.parse_args()

    if not isinstance(args.model, list):
        args.model = [args.model]


    # Now add all to a simple plot
    fig, axes = plt.subplots(1, 1, sharex=True, figsize=(5,5), dpi=100, facecolor='w', edgecolor='k')
    ax = axes
    ax.set_xlabel(r'$M_{\rm vir, DM} [M_{\odot}]$', fontsize=16)
    ax.set_ylabel(r'$M_{*} [M_{\odot}]$', fontsize=16)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.tick_params(labelsize=12)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.minorticks_on()
    ax.tick_params(which='both',axis="both",direction="in")

    # Add model to plot

    mhalo = np.linspace(10.5,13,100)
    mhalo = 10**mhalo
    if args.ref == 'Moster+2013':
        m_model = Moster2013(args.z,mhalo)

    ax.plot(mhalo,m_model, 'k-')

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
        os.chdir(simfolder)
        result = subprocess.run('IDtoZetas.out -ask '+str(int(args.z)), shell=True,stdout=subprocess.PIPE)
        os.chdir('../')
        indexout = int(result.stdout.decode('utf-8').split('is')[1].split('(z')[0])

        ozyfile = 'ozy_%05d.hdf5' % (indexout)
        sim = ozy.load(os.path.join(groupspath, ozyfile))
        physics = {'hydro':False,
                'metals':False,
                'magnetic':False,
                'cr':False,
                'rt':False,
                'bh':False,
                'AGN':False,
                'dust':False}
        varIDs = io_ramses.hydroID()
        io_ramses.read_hydrofile_descriptor(sim.simulation.fullpath,varIDs)
        if varIDs.density != 0 and varIDs.vx != 0:
            physics['hydro'] = True
        if varIDs.metallicity != 0:
            physics['metals'] = True
        if varIDs.blx != 0:
            physics['magnetic'] = True
        if varIDs.ecr != 0:
            physics['cr'] = True
        if varIDs.xhii != 0:
            physics['rt'] = True

        progind = args.ind
        if args.NUT:
            virial_mass = [i.virial_quantities['mass'] for i in sim.galaxies]
            progind = np.argmax(virial_mass)
        
        mstellar = sim.galaxies[progind].mass['stellar']
        mstellar = YTQuantity(mstellar, 'code_mass', registry=sim.unit_registry)
        mhalo = sim.galaxies[progind].halo.virial_quantities['mass']
        mhalo = YTQuantity(mhalo, 'code_mass', registry=sim.unit_registry)
        if not physics['cr'] and not physics['magnetic']:
            ax.plot(mhalo.in_units('Msun').d,mstellar.in_units('Msun').d, marker='s', markeredgecolor='b', markerfacecolor='none')
        elif not physics['cr'] and physics['magnetic']:
            ax.plot(mhalo.in_units('Msun').d,mstellar.in_units('Msun').d, marker='s', markeredgecolor='m', markerfacecolor='none')
        elif physics['cr']:
            ax.plot(mhalo.in_units('Msun').d,mstellar.in_units('Msun').d, marker='s', markeredgecolor='g', markerfacecolor='none')
    fig.subplots_adjust(top=0.91, bottom=0.1,right=0.97,hspace=0.0)
    if args.NUT:
        args.ind = 'NUT'
    fig.savefig(os.getcwd()+'/smf_'+str(args.ref)+'_'+str(args.ind)+'.png', format='png', dpi=300)