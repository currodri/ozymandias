"""
This code allows the comparison of the measured Stellar-to-Halo Mass Function in a cosmological simulation
to the prediction of models.

By: F. Rodriguez Montero (28/10/2021)
"""

# Import required libraries
import argparse
import itertools
import os
import subprocess
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import ozy
import seaborn as sns
from amr2 import io_ramses
from astropy.cosmology import FlatLambdaCDM, z_at_value
from matplotlib.collections import LineCollection

from numpy.core.fromnumeric import var
from ozy.utils import RotationAwareAnnotation

sns.set(style="white")
plt.rcParams["axes.axisbelow"] = False


def Moster2013(z, mhalo):
    """
    This function returns the best fit results by Moster et al. (2013) for a given redshift z.

    Throughout this paper, we assume a 7-year Wilkinson Mi- crowave Anisotropy Probe (WMAP7) LambdaCDM cosmology with
    (omega_m , omega_lambda , omega_b , h, n, sigma_8 ) = (0.272, 0.728, 0.046, 0.704, 0.967, 0.810). 
    We employ a Chabrier (2003) initial mass function (IMF) and we convert all stellar masses to this IMF

    Masses are returned in units of Msun.
    """
    fb = 0.046 / 0.272

    def compute_smhm(params,z,mhalo):
        rf = z / (1+z)
        # Equation 11 in Moster et al. (2013)
        M1 = 10**(params['M10'] + params['M11']*rf)
        
        # Equation 12 in Moster et al. (2013)
        N = params['N10'] + params['N11'] * rf

        # Equation 13 in Moster et al. (2013)
        beta = params['B10'] + params['B11'] * rf

        # Equation 14 in Moster et al. (2013)
        gamma= params['G10'] + params['G11'] * rf

        # Equation 2 in Moster et al. (2013)
        M_star = 2 * N * ((mhalo / M1)**(-beta) + (mhalo / M1)**gamma)**(-1.0)
        M_star = M_star * mhalo
        return M_star

    # Taken from Table 1 of Moster et al. (2013)
    fit_params = dict(
        M10 = np.array([11.590,0.236]),
        M11 = np.array([1.195,0.353]),
        N10 = np.array([0.0351,0.0058]),
        N11 = np.array([-0.0247,0.0069]),
        B10 = np.array([1.376,0.153]),
        B11 = np.array([-0.826,0.225]),
        G10 = np.array([0.608,0.059]),
        G11 = np.array([0.329,0.173])
        )
    fit_combs = [(item[0]+item[1],item[0]-item[1]) for item in fit_params.values()]

    if not os.path.exists('moster2013_low.npy') or not os.path.exists('moster2013_high.npy'):
        # Now, iterate over all combinations of parameters
        M_star_low = np.full(len(mhalo),np.infty)
        M_star_high = np.zeros(len(mhalo))
        ncomb = 0
        for comb in itertools.product(*fit_combs):
            p = {'M10':comb[0],
                'M11':comb[1],
                'N10':comb[2],
                'N11':comb[3],
                'B10':comb[4],
                'B11':comb[5],
                'G10':comb[6],
                'G11':comb[7],
            }
            M_star = compute_smhm(p,z,mhalo)
            for m in range(0, len(mhalo)):
                M_star_low[m] = min(M_star_low[m],M_star[m])
                M_star_high[m] = max(M_star_high[m],M_star[m])
            ncomb += 1
        print('ncomb: ',ncomb)
    
        np.save('moster2013_low.npy',M_star_low)
        np.save('moster2013_high.npy',M_star_high)
    else:
        print('Found Moster data!')
        M_star_low = np.load('moster2013_low.npy')
        M_star_high = np.load('moster2013_high.npy')

    # Label
    label = 'Moster et al. 2013 \n (z=%.2f)'%z
    return M_star_low,M_star_high,label

def Behroozi2013(z,mhalo):
    """
    This function returns the best fit results by Behroozi et al. (2013) for a given redshift z.

    We additionally assume a flat, ΛCDM cosmology with parameters ΩM = 0.27, ΩΛ = 0.73, h = 0.7, ns = 0.95, and sigma_8 = 0.82

    Masses are returned in units of Msun.
    """

    fb = 0.046 / 0.27

    def compute_smhm(params,a,z,mhalo):

        # Equation 4 (they're all together)
        nu = np.exp(-4.0*a**2)
        M1 = params['M_0'] + (params['M_a']*(a-1.) + params['M_z']*z)*nu
        M1 = 10**M1
        epsilon = params['epsilon_0'] + (params['epsilon_a']*(a-1.) + params['epsilon_z']*z)*nu + params['epsilon_a2']*(a-1.)
        epsilon = 10**epsilon
        alpha = params['alpha_0'] + params['alpha_a']*(a-1.)*nu
        delta = params['delta_0'] + (params['delta_a']*(a-1.) + params['delta_z']*z)*nu
        gamma = params['gamma_0'] + (params['gamma_a']*(a-1.) + params['gamma_z']*z)*nu
        # Final equations
        x = np.log10(mhalo/M1)
        f0 = -np.log10(2.0) + delta * (np.log10(2.)**gamma/(1+np.e))
        f = -np.log10(10**(alpha*x)+1) + delta*(np.log10(1+np.exp(x)))**gamma/(1.+np.exp(10**(-x)))
        M_star = np.log10(epsilon*M1) + f - f0
        M_star = 10**M_star

        return M_star

    # Taken from the first fit in Table J1 (i.e. Obs, all galaxies, central and satellites)
    # First column is the parameters for the best-fitting model, and the second and third are the
    # 68% confidence interval for the model posterior distribution
    fit_params = dict(
        epsilon_0 = np.array([-1.777,0.133,-0.146]),
        alpha_0 = np.array([-1.412,0.020,-0.105]),
        epsilon_a = np.array([-0.006,0.113,-0.361]),
        alpha_a   = np.array([0.731,0.344,-0.296]),
        epsilon_a2 = np.array([-0.119,0.061,-0.012]),
        epsilon_z = np.array([0.000,0.003,-0.104]),
        M_0 = np.array([11.514,0.053,-0.009]),
        M_a = np.array([-1.793,0.315,-0.330]),
        delta_0 = np.array([3.508,0.087,-0.369]),
        delta_a = np.array([2.608,2.446,-1.261]),
        delta_z = np.array([-0.043,0.958,-0.071]),
        gamma_0 = np.array([0.316,0.076,-0.012]),
        M_z = np.array([-0.731,0.464,-0.064]),
        gamma_a = np.array([1.319,0.584,-0.505]),
        gamma_z = np.array([0.279,0.256,-0.081]),
    )
    fit_combs = [(item[0]+item[1],item[0]+item[2]) for item in fit_params.values()]
    a = 1. / (1. + z)
    if not os.path.exists('behroozi2013_low.npy') or not os.path.exists('behroozi2013_high.npy'):
        # Now, iterate over all combinations of parameters
        M_star_low = np.full(len(mhalo),np.infty)
        M_star_high = np.zeros(len(mhalo))
        ncomb = 0
        for comb in itertools.product(*fit_combs):
            p = {'epsilon_0':comb[0],
                'alpha_0':comb[1],
                'epsilon_a':comb[2],
                'alpha_a':comb[3],
                'epsilon_a2':comb[4],
                'epsilon_z':comb[5],
                'M_0':comb[6],
                'M_a':comb[7],
                'delta_0':comb[8],
                'delta_a':comb[9],
                'delta_z':comb[10],
                'gamma_0':comb[11],
                'M_z':comb[12],
                'gamma_a':comb[13],
                'gamma_z':comb[14],
            }
            M_star = compute_smhm(p,a,z,mhalo)
            for m in range(0, len(mhalo)):
                M_star_low[m] = min(M_star_low[m],M_star[m])
                M_star_high[m] = max(M_star_high[m],M_star[m])
            ncomb += 1
        print('ncomb: ',ncomb)
    
        np.save('behroozi2013_low.npy',M_star_low)
        np.save('behroozi2013_high.npy',M_star_high)
    else:
        print('Found Behroozi data!')
        M_star_low = np.load('behroozi2013_low.npy')
        M_star_high = np.load('behroozi2013_high.npy')

    label = 'Behroozi et al. 2013 (z=%.2f)'%z

    return M_star_low,M_star_high,label

def Behroozi2019(z,mhalo):
    """
    This function returns the best fit results by Behroozi et al. (2019) for a given redshift z.

    We adopt a flat, LambdaCDM cosmology with parameters (omega_m = 0.307, omega_lambda = 0.693, h = 0.678, sigma_8 = 0.823, ns = 0.96) 
    consistent with Planck results (Planck Collaboration et al. 2016).

    Masses are returned in units of Msun.
    """

    fb = 0.16061638603207626

    def compute_smhm(params,a,z,mhalo):
        # Equation J3
        M1 = params['M_0'] + params['M_a']*(a-1.) - params['M_lna']*np.log(a) + params['M_z']*z 
        M1 = 10**M1
        # Equation J4
        epsilon = params['epsilon_0'] + params['epsilon_a']*(a-1.) - params['epsilon_lna']*np.log(a) + params['epsilon_z']*z
        # Equation J5
        alpha = params['alpha_0'] + params['alpha_a']*(a-1.) - params['alpha_lna']*np.log(a) + params['alpha_z']*z
        # Equation J6
        beta = params['beta_0'] + params['beta_a']*(a-1.) + params['beta_z']*z
        # Equation J7
        delta = params['delta_0']
        # Equation J8
        gamma = params['gamma_0'] + params['gamma_a']*(a-1.) + params['gamma_z']*z
        gamma = 10**gamma
        # Final equations
        # Equation J2
        x = np.log10(mhalo/M1)
        # Equation J1
        M_star = epsilon - np.log10(10**(-alpha*x) + 10**(-beta*x)) + gamma*np.exp(-0.5*(x/delta)**2)
        M_star = M1 * 10**M_star

        return M_star

    # Taken from the first fit in Table J1 (i.e. Obs, all galaxies, central and satellites)
    # First column is the parameters for the best-fitting model, and the second and third are the
    # 68% confidence interval for the model posterior distribution
    fit_params = dict(
        epsilon_0 = np.array([-1.435,0.023,-0.075]),
        alpha_lna = np.array([-1.732,0.519,-1.271]),
        epsilon_a = np.array([1.831,-0.066,-2.818]),
        alpha_z   = np.array([0.178,0.196,-0.084]),
        epsilon_lna = np.array([1.368,0.089,-2.557]),
        beta_0 = np.array([0.482,0.053,0.002]),
        epsilon_z = np.array([-0.217,0.466,-0.048]),
        beta_a = np.array([-0.841,0.700,0.188]),
        M_0 = np.array([12.035,0.008,-0.100]),
        beta_z = np.array([-0.471,0.288,0.049]),
        M_a = np.array([4.556,0.076,-2.453]),
        delta_0 = np.array([0.411,0.060,-0.087]),
        M_lna = np.array([4.417,0.237,-2.255]),
        gamma_0 = np.array([-1.034,0.264,-0.177]),
        M_z = np.array([-0.731,0.464,-0.064]),
        gamma_a = np.array([-3.100,2.496,-0.363]),
        alpha_0 = np.array([1.963,0.163,-0.010]),
        gamma_z = np.array([-1.055,0.973,-0.127]),
        alpha_a = np.array([-2.316,0.855,-1.361]),
    )
    fit_combs = [(item[0]+item[1],item[0]+item[2]) for item in fit_params.values()]
    a = 1. / (1. + z)
    if not os.path.exists('behroozi2019_low.npy') or not os.path.exists('behroozi2019_high.npy'):
        # Now, iterate over all combinations of parameters
        M_star_low = np.full(len(mhalo),np.infty)
        M_star_high = np.zeros(len(mhalo))
        ncomb = 0
        for comb in itertools.product(*fit_combs):
            p = {'epsilon_0':comb[0],
                'alpha_lna':comb[1],
                'epsilon_a':comb[2],
                'alpha_z':comb[3],
                'epsilon_lna':comb[4],
                'beta_0':comb[5],
                'epsilon_z':comb[6],
                'beta_a':comb[7],
                'M_0':comb[8],
                'beta_z':comb[9],
                'M_a':comb[10],
                'delta_0':comb[11],
                'M_lna':comb[12],
                'gamma_0':comb[13],
                'M_z':comb[14],
                'gamma_a':comb[15],
                'alpha_0':comb[16],
                'gamma_z':comb[17],
                'alpha_a':comb[18]
            }
            M_star = compute_smhm(p,a,z,mhalo)
            for m in range(0, len(mhalo)):
                M_star_low[m] = min(M_star_low[m],M_star[m])
                M_star_high[m] = max(M_star_high[m],M_star[m])
            ncomb += 1
        print('ncomb: ',ncomb)
    
        np.save('behroozi2019_low.npy',M_star_low)
        np.save('behroozi2019_high.npy',M_star_high)
    else:
        print('Found Behroozi data!')
        M_star_low = np.load('behroozi2019_low.npy')
        M_star_high = np.load('behroozi2019_high.npy')

    label = 'Behroozi et al. 2019 (z=%.2f)'%z

    return M_star_low,M_star_high,label


if __name__ == '__main__':

    # Parse the command line arguments.
    parser = argparse.ArgumentParser(description='Tracking of galaxy masses across cosmic time.')
    parser.add_argument('model', type=str, nargs='+', help='Simulation names from which extract OZY results.')
    parser.add_argument('--ind', type=int, default=0, help='Index of the galaxy in the ending snapshot.')
    parser.add_argument('--z',type=float,default=[3.0],nargs='+', help="Redshift at which to extract data.")
    parser.add_argument('--ref', type=str, default='all', help='Model to which compare.')
    parser.add_argument('--NUT', type=bool, default=True, help='If True, it looks for NUT as the most massive galaxy (stars) in the last snapshot.')
    args = parser.parse_args()

    if not isinstance(args.model, list):
        args.model = [args.model]


    # Now add all to a simple plot
    fig, axes = plt.subplots(1, 1,figsize=(7,7), dpi=100, facecolor='w', edgecolor='k')
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
    minmass = np.log10(3e+9)
    maxmass = np.log10(4e+11)
    ax.set_xlim([10**minmass,10**maxmass])
    ax.set_ylim([1e6,1e+11])
    mhalo = np.linspace(minmass,maxmass,100)
    mhalo = 10**mhalo
    if args.ref == 'Moster+2013':
        m_model_low,m_model_high,m_label = Moster2013(0.0,mhalo)

        ax.fill_between(x=mhalo,y1=m_model_low,y2=m_model_high,color='b',alpha=0.2)
        ax.plot(mhalo,m_model_low,color='b')
        ax.plot(mhalo,m_model_high,color='b')
        ind = mhalo.searchsorted(1.5e+10)
        ax.text(1.5e+10,0.25*(m_model_high[ind]+m_model_low[ind]), m_label, fontsize=16, color='b',
                rotation=45,horizontalalignment='center', verticalalignment='center')
    elif args.ref == 'Behroozi+2019':
        m_model_low,m_model_high,m_label = Behroozi2019(0.0,mhalo)

        ax.fill_between(x=mhalo,y1=m_model_low,y2=m_model_high,color='orange',alpha=0.2)
        ax.plot(mhalo,m_model_low,color='orange')
        ax.plot(mhalo,m_model_high,color='orange')
        ind = mhalo.searchsorted(8e+10)
        ax.text(8e+10,0.4*(m_model_high[ind]+m_model_low[ind]), m_label, fontsize=16, color='orange',
                rotation=40,horizontalalignment='center', verticalalignment='center')
    elif args.ref == 'Behroozi+2013':
        m_model_low,m_model_high,m_label = Behroozi2013(0.0,mhalo)

        ax.fill_between(x=mhalo,y1=m_model_low,y2=m_model_high,color='orange',alpha=0.2)
        ax.plot(mhalo,m_model_low,color='orange')
        ax.plot(mhalo,m_model_high,color='orange')
        ind = mhalo.searchsorted(8e+10)
        ax.text(8e+10,0.4*(m_model_high[ind]+m_model_low[ind]), m_label, fontsize=16, color='orange',
                rotation=40,horizontalalignment='center', verticalalignment='center')
    elif args.ref == 'all':
        m_model_low,m_model_high,m_label = Moster2013(0.0,mhalo)

        ax.fill_between(x=mhalo,y1=m_model_low,y2=m_model_high,color='b',alpha=0.2)
        ax.plot(mhalo,m_model_low,color='b')
        ax.plot(mhalo,m_model_high,color='b')
        ind = mhalo.searchsorted(6e+10)
        ax.text(6e+10,0.25*(m_model_high[ind]+m_model_low[ind]), m_label, fontsize=12, color='b',
                rotation=45,horizontalalignment='center', verticalalignment='center')
            
        m_model_low,m_model_high,m_label = Behroozi2013(0.0,mhalo)

        ax.fill_between(x=mhalo,y1=m_model_low,y2=m_model_high,color='orange',alpha=0.2)
        ax.plot(mhalo,m_model_low,color='orange')
        ax.plot(mhalo,m_model_high,color='orange')
        ind = mhalo.searchsorted(6e+9)
        ax.text(8e+9,0.5*(m_model_high[ind]+m_model_low[ind]), m_label, fontsize=12, color='orange',
                rotation=30,horizontalalignment='center', verticalalignment='center')
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
        mstellar = np.zeros(len(args.z))
        mhalo = np.zeros(len(args.z))
        for j in range(0, len(args.z)):
            groupspath = os.path.join(simfolder, 'Groups')
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
            
            mstellar[j] = sim.galaxies[progind].mass['stellar'].to('Msun')
            mhalo[j] = sim.galaxies[progind].halo.virial_quantities['mass'].to('Msun')
        if not sim.simulation.physics['cr'] and not sim.simulation.physics['magnetic']:
            ax.plot(mhalo,mstellar, markersize=10, marker="d", markeredgecolor='b',
                    markerfacecolor='b',label=args.model[i],color='b')
        elif not sim.simulation.physics['cr'] and sim.simulation.physics['magnetic']:
            ax.plot(mhalo,mstellar, markersize=10,marker='*', markeredgecolor='m',
                    markerfacecolor='m',label=args.model[i],color='m')
        elif sim.simulation.physics['cr']:
            ax.plot(mhalo,mstellar, markersize=10,marker='v', markeredgecolor='g',
                    markerfacecolor='g',label=args.model[i],color='g')
        if i == 0:
            # Add fixed fb to plot
            omega_dm = 1.0 - (sim.simulation.omega_lambda + sim.simulation.omega_k + sim.simulation.omega_baryon)
            fb = sim.simulation.omega_baryon / omega_dm
            mhalo_model = np.linspace(minmass,maxmass,100)
            mhalo_model = 10**mhalo_model
            mstar = lambda x: 0.3*x
            ax.plot(mhalo_model,mstar(mhalo_model),linestyle='dashdot',color='k')
            midx = 0.5*(maxmass+minmass)
            an = RotationAwareAnnotation(r'$M_*=f_{\rm b}M_{\rm DM}$', xy=(10**midx,mstar(1.1*10**midx)), p=(10**(midx+0.5),mstar(1.1*10**(midx+0.5))), ax=ax,
                                        xytext=(-1,1), textcoords="offset points", 
                                        ha="center", va="baseline", fontsize=16)
        # Compute conversion fractions just to print to screen (may be useful)
        print(20*'-')
        print('BARYON CONVERSION EFFICIENCY [%]')
        print('Measured for model %s at the redshift bins of '%(args.model[i]),args.z)
        max_mstellar = mstar(mhalo)
        efficiency = 100 * mstellar / max_mstellar
        print(efficiency)
        print(20*'-')

    # Plot vertical lines separating redshifts
    ax.text(4.7e+9,2e+6, 'z=8', fontsize=16, color='k',alpha=0.5)
    ax.plot([1e+10,1e+10],[1e6,1e+11],'k--',alpha=0.3) # z=6
    ax.text(2e+10,2e+6, 'z=6', fontsize=16, color='k',alpha=0.5)
    ax.plot([5e+10,5e+10],[1e6,1e+11],'k--',alpha=0.3) # z=3
    ax.text(6e+10,2e+6, 'z=3', fontsize=16, color='k',alpha=0.5)
    # ax.plot([4e+9,1e+10],[1e6,1e+11],'k--',alpha=0.3) # z=3.5
    ax.plot([1.2e+11,1.2e+11],[1e6,1e+11],'k--',alpha=0.3) # z=2
    ax.text(1.7e+11,2e+6, 'z=2', fontsize=16, color='k',alpha=0.5)

    fig.subplots_adjust(top=0.91, bottom=0.1,right=0.97,hspace=0.0)
    if args.NUT:
        args.ind = 'NUT'
    fig.savefig(os.getcwd()+'/smf_'+str(args.ref)+'_'+str(args.ind)+'.png', format='png', dpi=300)
