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
from uncertainties import ufloat
from uncertainties.umath import *

sns.set(style="white")
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "Computer Modern Roman",
})
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
        M10 = ufloat(params['M10'][0],params['M10'][1])
        M11 = ufloat(params['M11'][0],params['M11'][1])
        M1 = M10 + rf*M11
        M1 = pow(10,M1)
        
        # Equation 12 in Moster et al. (2013)
        N10 = ufloat(params['N10'][0],params['N10'][1])
        N11 = ufloat(params['N11'][0],params['N11'][1])
        N = N10 + N11 * rf

        # Equation 13 in Moster et al. (2013)
        B10 = ufloat(params['B10'][0],params['B10'][1])
        B11 = ufloat(params['B11'][0],params['B11'][1])
        beta = B10 + B11 * rf

        # Equation 14 in Moster et al. (2013)
        G10 = ufloat(params['G10'][0],params['G10'][1])
        G11 = ufloat(params['G11'][0],params['G11'][1])
        gamma = G10 + G11 * rf

        # Equation 2 in Moster et al. (2013)
        M_star = 2 * N * pow(pow(mhalo / M1,-beta) + pow(mhalo / M1,gamma),-1)
        M_star = M_star * mhalo

        return M_star.n,M_star.s

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

    # Compute SMHM
    M_star,delta_M_star = np.zeros(len(mhalo)),np.zeros(len(mhalo))
    for m in range(0,len(mhalo)):
        M_star[m],delta_M_star[m] = compute_smhm(fit_params,z,mhalo[m])
    # Label
    label = 'Moster et al. 2013 \n'+ r'($z=%.2f$)'%z
    return M_star-delta_M_star,M_star+delta_M_star,label

# def Moster2018(z,mhalo):
#     """
#     This function returns the best fit results by Moster et al. (2018) for a given redshift z.

#     Masses are returned in units of Msun.
#     """
    
#     # Fitting results from MCMC (Table 6)
#     fit_params = dict(
#         M0 = np.array([11.339,0.005,0.008]),
#         MZ = np.array([0.692,0.010,0.009]),
#         eps0 = np.array([0.005,0.001,0.001]),
#         epsZ = np.array([0.689,0.003,0.003]),
#         beta0 = np.array([3.344,0.084,0.101]),
#         betaZ = np.array([-2.079,0.127,0.134]),
#         gamma0 = np.array([0.600,0.002,0.003])
#         )
    
#     def compute_smhm(params,z,mhalo,err=0):
#         rf = z / (1+z)
#         # Equation 5 in Moster et al. (2018)
#         M0 = ufloat(params['M0'][0],params['M0'][1+err])
#         MZ = ufloat(params['MZ'][0],params['MZ'][1+err])
#         M1 = M0 + rf*MZ
#         M1 = pow(10,M1)
        
#         # Equation 7 in Moster et al. (2018)
#         eps0 = ufloat(params['eps0'][0],params['eps0'][1+err])
#         epsZ = ufloat(params['epsZ'][0],params['epsZ'][1+err])
#         epsN = eps0 + epsZ * rf

#         # Equation 8 in Moster et al. (2018)
#         beta0 = ufloat(params['beta0'][0],params['beta0'][1+err])
#         betaZ = ufloat(params['betaZ'][0],params['betaZ'][1+err])
#         beta = beta0 + betaZ * rf

#         # Equation 10 in Moster et al. (20138)
#         gamma = ufloat(params['gamma0'][0],params['gamma0'][1+err])
#         print(gamma)
#         # Equation 2 in Moster et al. (2013)
#         M_star = 2 * epsN * pow(pow(mhalo / M1,-beta) + pow(mhalo / M1,gamma),-1)
#         M_star = M_star * mhalo

#         return M_star.n,M_star.s
    
#     # Compute SMHM
#     M_star,delta_plus,delta_minus = np.zeros(len(mhalo)),np.zeros(len(mhalo)),np.zeros(len(mhalo))
#     for m in range(0,len(mhalo)):
#         M_star[m],delta_plus[m] = compute_smhm(fit_params,z,mhalo[m])
#         M_star[m],delta_minus[m] = compute_smhm(fit_params,z,mhalo[m],err=1)
#     # Label
#     label = 'Moster et al. 2018 \n'+ r'($z=%.2f$)'%z
#     return M_star-delta_minus,M_star+delta_plus,label
    
def Moster2018(z,mhalo):
    """
    This function returns the best fit results by Moster et al. (2018) for a given redshift z.
    Masses are returned in units of Msun.
    """
    # Fitting parameters for all galaxies (Table 8)
    
    fb = 0.156
    
    redshifts = np.array([0.1,0.5,1.0,2.0,4.0,8.0])
    fit_params = {
        'M1':[11.78,11.86,11.98,11.99,12.07,12.10],
        'epsN':[0.15,0.18,0.19,0.19,0.20,0.24],
        'beta':[1.78,1.67,1.53,1.46,1.36,1.30],
        'gamma':[0.57,0.58,0.59,0.59,0.59,0.60,0.60],
        'Msigma':[10.85,10.80,10.75,10.70,10.60,10.40],
        'sigma0':[0.16,0.14,0.12,0.10,0.06,0.02],
        'alpha':[1.00,0.75,0.60,0.45,0.35,0.30]
    }
    pos = np.searchsorted(redshifts,z)
    def efficiency(mh):
        eff = 2. * fit_params['epsN'][pos]*(1./((mh/10**fit_params['M1'][pos])**(-fit_params['beta'][pos]) + (mh/10**fit_params['M1'][pos])**(fit_params['gamma'][pos])))
        return eff
    def scatter(mh):
        sigma = fit_params['sigma0'][pos] + np.log10((mh/10**fit_params['Msigma'][pos])**(-fit_params['alpha'][pos])+1.)
        return sigma
    M_star,delta_M_star = np.zeros(len(mhalo)),np.zeros(len(mhalo))
    for m in range(0,len(mhalo)):
        
        M_star[m] = mhalo[m] * efficiency(mhalo[m])
        delta_M_star[m] = scatter(mhalo[m])
    # Label
    label = 'Moster et al. 2018 \n'+ r'($z=%.2f$)'%redshifts[pos]
    return 10**(np.log10(M_star)-delta_M_star),10**(np.log10(M_star)+delta_M_star),label


# def Behroozi2013(z,mhalo):
#     """
#     This function returns the best fit results by Behroozi et al. (2013) for a given redshift z.

#     We additionally assume a flat, ΛCDM cosmology with parameters ΩM = 0.27, ΩΛ = 0.73, h = 0.7, ns = 0.95, and sigma_8 = 0.82

#     Masses are returned in units of Msun.
#     """

#     fb = 0.046 / 0.27

#     def mass(params,a,z,mhalo):

#         # Equation 4 (they're all together)
#         nu = np.exp(-4.0*a**2)
#         M1 = params['M_0'][0] + (params['M_a'][0]*(a-1.) + params['M_z'][0]*z)*nu
#         M1 = 10**M1
#         epsilon = params['epsilon_0'][0] + (params['epsilon_a'][0]*(a-1.) + params['epsilon_z'][0]*z)*nu + params['epsilon_a2'][0]*(a-1.)
#         epsilon = 10**epsilon
#         alpha = params['alpha_0'][0] + params['alpha_a'][0]*(a-1.)*nu
#         delta = params['delta_0'][0] + (params['delta_a'][0]*(a-1.) + params['delta_z'][0]*z)*nu
#         gamma = params['gamma_0'][0] + (params['gamma_a'][0]*(a-1.) + params['gamma_z'][0]*z)*nu
#         # Final equations
#         x = np.log10(mhalo/M1)
#         f0 = -np.log10(2.0) + delta * (np.log10(2.)**gamma/(1+np.e))
#         f = -np.log10(10**(alpha*x)+1) + delta*(np.log10(1+np.exp(x)))**gamma/(1.+np.exp(10**(-x)))
#         M_star = np.log10(epsilon*M1) + f - f0
#         M_star = 10**M_star

#         return M_star
    
#     def scatter(params,a,z,mhalo):
#         # Equation 11
#         eta = params['eta_0'][0] + params['eta_z'][0]*(a-0.1)
#         return eta

#     # Taken from the equations in page 9 and 10
#     fit_params = dict(
#         epsilon_0 = np.array([-1.777,0.133,-0.146]),
#         alpha_0 = np.array([-1.412,0.020,-0.105]),
#         epsilon_a = np.array([-0.006,0.113,-0.361]),
#         alpha_a   = np.array([0.731,0.344,-0.296]),
#         epsilon_a2 = np.array([-0.119,0.061,-0.012]),
#         epsilon_z = np.array([0.000,0.003,-0.104]),
#         M_0 = np.array([11.514,0.053,-0.009]),
#         M_a = np.array([-1.793,0.315,-0.330]),
#         delta_0 = np.array([3.508,0.087,-0.369]),
#         delta_a = np.array([2.608,2.446,-1.261]),
#         delta_z = np.array([-0.043,0.958,-0.071]),
#         gamma_0 = np.array([0.316,0.076,-0.012]),
#         M_z = np.array([-0.731,0.464,-0.064]),
#         gamma_a = np.array([1.319,0.584,-0.505]),
#         gamma_z = np.array([0.279,0.256,-0.081]),
#         eta_0 = np.array([0.218,0.011,0.033]),
#         eta_z = np.array([-0.023,0.052,0.068])
#     )
#     a = 1. / (1. + z)
#     M_star,delta_M_star = np.zeros(len(mhalo)),np.zeros(len(mhalo))
#     for m in range(0,len(mhalo)):
#         M_star[m] = mass(fit_params,a,z,mhalo[m])
#         delta_M_star[m] = scatter(fit_params,a,z,mhalo[m])
#     label = 'Behroozi et al. 2013'+ r'($z=%.2f$)'%z

#     return 10**(np.log10(M_star)-delta_M_star),10**(np.log10(M_star)+delta_M_star),label
def Behroozi2013(z,mhalo):
    """
    This function returns the best fit results by Behroozi et al. (2013) for a given redshift z.

    We additionally assume a flat, ΛCDM cosmology with parameters ΩM = 0.27, ΩΛ = 0.73, h = 0.7, ns = 0.95, and sigma_8 = 0.82

    Masses are returned in units of Msun.
    """

    fb = 0.046 / 0.27

    def compute_smhm(params,a,z,mhalo,err=0):

        # Equation 4 (they're all together)
        nu = np.exp(-4.0*a**2)
        M0 = ufloat(params['M_0'][0],params['M_0'][1+err])
        Ma = ufloat(params['M_a'][0],params['M_a'][1+err])
        Mz = ufloat(params['M_z'][0],params['M_z'][1+err])
        M1 = M0 + (Ma*(a-1.) + Mz*z)*nu
        M1 = pow(10,M1)
        epsilon0 = ufloat(params['epsilon_0'][0],params['epsilon_0'][1+err])
        epsilona = ufloat(params['epsilon_a'][0],params['epsilon_a'][1+err])
        epsilonz = ufloat(params['epsilon_z'][0],params['epsilon_z'][1+err])
        epsilona2 = ufloat(params['epsilon_a2'][0],params['epsilon_a2'][1+err])

        epsilon = epsilon0 + (epsilona*(a-1.) + epsilonz*z)*nu + epsilona2*(a-1.)
        epsilon = pow(10,epsilon)
        alpha0 = ufloat(params['alpha_0'][0],params['alpha_0'][1+err])
        alphaa = ufloat(params['alpha_a'][0],params['alpha_a'][1+err])
        alpha = alpha0 + alphaa*(a-1.)*nu
        
        delta0 = ufloat(params['delta_0'][0],params['delta_0'][1+err])
        deltaa = ufloat(params['delta_a'][0],params['delta_a'][1+err])
        deltaz = ufloat(params['delta_z'][0],params['delta_z'][1+err])
        delta = delta0 + (deltaa*(a-1.) + deltaz*z)*nu
        
        gamma0 = ufloat(params['gamma_0'][0],params['gamma_0'][1+err])
        gammaa = ufloat(params['gamma_a'][0],params['gamma_a'][1+err])
        gammaz = ufloat(params['gamma_z'][0],params['gamma_z'][1+err])
        gamma = gamma0 + (gammaa*(a-1.) + gammaz*z)*nu
        # Final equations
        x = log10(mhalo/M1)
        f0 = -log10(2.0) + delta * (log10(2.)**gamma/(1+np.e))
        f = -log10(pow(10,alpha*x)+1) + delta*(log10(1+exp(x)))**gamma/(1.+exp(pow(10,-x)))
        M_star = log10(epsilon*M1) + f - f0
        M_star = pow(10,M_star)

        return M_star.n,M_star.s
    
    def scatter(params,a,z,mhalo):
        # Equation 11
        eta = params['eta_0'][0] + params['eta_z'][0]*(a-0.1)
        return eta

    # Taken from the equations in page 9 and 10
    fit_params = dict(
        epsilon_0 = np.array([-1.777,0.133,0.146]),
        alpha_0 = np.array([-1.412,0.020,0.105]),
        epsilon_a = np.array([-0.006,0.113,0.361]),
        alpha_a   = np.array([0.731,0.344,0.296]),
        epsilon_a2 = np.array([-0.119,0.061,0.012]),
        epsilon_z = np.array([0.000,0.003,0.104]),
        M_0 = np.array([11.514,0.053,0.009]),
        M_a = np.array([-1.793,0.315,0.330]),
        delta_0 = np.array([3.508,0.087,0.369]),
        delta_a = np.array([2.608,2.446,1.261]),
        delta_z = np.array([-0.043,0.958,0.071]),
        gamma_0 = np.array([0.316,0.076,0.012]),
        M_z = np.array([-0.731,0.464,0.064]),
        gamma_a = np.array([1.319,0.584,0.505]),
        gamma_z = np.array([0.279,0.256,0.081]),
        eta_0 = np.array([0.218,0.011,0.033]),
        eta_z = np.array([-0.023,0.052,0.068])
    )
    # Compute SMHM
    M_star,delta_plus,delta_minus,delta_M_star = np.zeros(len(mhalo)),np.zeros(len(mhalo)),np.zeros(len(mhalo)),np.zeros(len(mhalo))
    a = 1. / (1. + z)
    for m in range(0,len(mhalo)):
        M_star[m],delta_plus[m] = compute_smhm(fit_params,a,z,mhalo[m])
        M_star[m],delta_minus[m] = compute_smhm(fit_params,a,z,mhalo[m],err=1)
        delta_M_star[m] = scatter(fit_params,a,z,mhalo[m])

    label = r'Behroozi et al. 2013 ($z=%.2f$)'%z
    print(delta_minus,delta_plus)
    M_down = M_star-delta_minus
    M_up = M_star+delta_plus
    return 10**(np.log10(M_down)-delta_M_star),10**(np.log10(M_up)+delta_M_star),label

def Behroozi2019(z,mhalo):
    """
    This function returns the best fit results by Behroozi et al. (2019) for a given redshift z.

    We adopt a flat, LambdaCDM cosmology with parameters (omega_m = 0.307, omega_lambda = 0.693, h = 0.678, sigma_8 = 0.823, ns = 0.96) 
    consistent with Planck results (Planck Collaboration et al. 2016).

    Masses are returned in units of Msun.
    """

    fb = 0.16061638603207626

    def compute_smhm(params,a,z,mhalo,err=0):
        # Equation J3
        M0 = ufloat(params['M_0'][0], params['M_0'][1+err])
        Ma = ufloat(params['M_a'][0], params['M_a'][1+err])
        Mlna = ufloat(params['M_lna'][0], params['M_lna'][1+err])
        Mz = ufloat(params['M_z'][0], params['M_z'][1+err])
        
        M1 = M0 + Ma*(a-1.) - Mlna*np.log(a) + Mz*z 
        M1 = pow(10,M1)
        
        # Equation J4
        eps0 = ufloat(params['epsilon_0'][0], params['epsilon_0'][1+err])
        epsa = ufloat(params['epsilon_a'][0], params['epsilon_a'][1+err])
        epslna = ufloat(params['epsilon_lna'][0], params['epsilon_lna'][1+err])
        epsz = ufloat(params['epsilon_z'][0], params['epsilon_z'][1+err])

        epsilon = eps0 + epsa*(a-1.) - epslna*np.log(a) + epsz*z
        
        # Equation J5
        alpha0 = ufloat(params['alpha_0'][0], params['alpha_0'][1+err])
        alphaa = ufloat(params['alpha_a'][0], params['alpha_a'][1+err])
        alphalna = ufloat(params['alpha_lna'][0], params['alpha_lna'][1+err])
        alphaz = ufloat(params['alpha_z'][0], params['alpha_z'][1+err])
        
        alpha = alpha0 + alphaa*(a-1.) - alphalna*np.log(a) + alphaz*z
        
        # Equation J6
        beta0 = ufloat(params['beta_0'][0], params['beta_0'][1+err])
        betaa = ufloat(params['beta_a'][0], params['beta_a'][1+err])
        betaz = ufloat(params['beta_z'][0], params['beta_z'][1+err])
        
        beta = beta0 + betaa*(a-1.) + betaz*z
        
        # Equation J7
        delta = ufloat(params['delta_0'][0],params['delta_0'][1+err])
        
        # Equation J8
        gamma0 = ufloat(params['gamma_0'][0], params['gamma_0'][1+err])
        gammaa = ufloat(params['gamma_a'][0], params['gamma_a'][1+err])
        gammaz = ufloat(params['gamma_z'][0], params['gamma_z'][1+err])
        gamma = gamma0 + gammaa*(a-1.) + gammaz*z
        gamma = pow(10,gamma)
        print(M1,epsilon,alpha,beta,delta,gamma)
        # Final equations
        # Equation J2
        x = log10(mhalo/M1)
        # Equation J1
        M_star = epsilon - log10(pow(10,-alpha*x) + pow(10,-beta*x)) + gamma*exp(-0.5*(x/delta)**2)
        M_star = M1 * pow(10,M_star)
        return M_star.n,M_star.s

    # Taken from the first fit in Table J1 (i.e. Obs, all galaxies, central and satellites)
    # First column is the parameters for the best-fitting model, and the second and third are the
    # 68% confidence interval for the model posterior distribution
    # fit_params = dict(
    #     epsilon_0 = np.array([-1.435,0.023,0.075]),
    #     alpha_lna = np.array([-1.732,0.519,1.271]),
    #     epsilon_a = np.array([1.831,0.066,0.066]),
    #     alpha_z   = np.array([0.178,0.084,0.084]),
    #     epsilon_lna = np.array([1.368,0.089,0.089]),
    #     beta_0 = np.array([0.482,0.053,0.002]),
    #     epsilon_z = np.array([-0.217,0.048,0.048]),
    #     beta_a = np.array([-0.841,0.188,0.188]),
    #     M_0 = np.array([12.035,0.008,0.008]),
    #     beta_z = np.array([-0.471,0.049,0.049]),
    #     M_a = np.array([4.556,0.076,0.076]),
    #     delta_0 = np.array([0.411,0.060,0.087]),
    #     M_lna = np.array([4.417,0.237,0.237]),
    #     gamma_0 = np.array([-1.034,0.264,0.177]),
    #     M_z = np.array([-0.731,0.064,0.064]),
    #     gamma_a = np.array([-3.100,0.363,0.363]),
    #     alpha_0 = np.array([1.963,0.163,0.010]),
    #     gamma_z = np.array([-1.055,0.127,0.127]),
    #     alpha_a = np.array([-2.316,0.855,1.361]),
    # )
    # fit_params = dict(
    #     epsilon_0 = np.array([-1.435,0.023,0.075]),
    #     alpha_lna = np.array([-1.732,0.519,1.271]),
    #     epsilon_a = np.array([1.831,0.066,2.818]),
    #     alpha_z   = np.array([0.178,0.196,0.084]),
    #     epsilon_lna = np.array([1.368,0.089,2.557]),
    #     beta_0 = np.array([0.482,0.053,0.002]),
    #     epsilon_z = np.array([-0.217,0.466,0.048]),
    #     beta_a = np.array([-0.841,0.700,0.188]),
    #     M_0 = np.array([12.035,0.008,0.100]),
    #     beta_z = np.array([-0.471,0.288,0.049]),
    #     M_a = np.array([4.556,0.076,2.453]),
    #     delta_0 = np.array([0.411,0.060,0.087]),
    #     M_lna = np.array([4.417,0.237,2.255]),
    #     gamma_0 = np.array([-1.034,0.264,0.177]),
    #     M_z = np.array([-0.731,0.464,0.064]),
    #     gamma_a = np.array([-3.100,2.496,0.363]),
    #     alpha_0 = np.array([1.963,0.163,0.010]),
    #     gamma_z = np.array([-1.055,0.973,0.127]),
    #     alpha_a = np.array([-2.316,0.855,1.361]),
    # )
    # fit_params = dict(
    #     epsilon_0 = np.array([-1.431495e+00,2.156496e-02,1.333165e-01]),
    #     alpha_lna = np.array([-1.816299e+00,4.326952e-01,1.388066e+00]),
    #     epsilon_a = np.array([1.757030e+00,-2.236675e-01,2.916873e+00]),
    #     alpha_z   = np.array([1.820800e-01,2.099754e-01,7.516579e-02]),
    #     epsilon_lna = np.array([1.350451e+00,-6.305647e-02,2.736450e+00]),
    #     beta_0 = np.array([4.702271e-01,6.133160e-02,-3.179759e-03]),
    #     epsilon_z = np.array([-2.178460e-01,4.804196e-01,3.468697e-02]),
    #     beta_a = np.array([-8.751643e-01,7.385737e-01,-2.069305e-01]),
    #     M_0 = np.array([1.207402e+01,1.219273e-02,8.281898e-02]),
    #     beta_z = np.array([-4.866420e-01,3.124524e-01,-6.955197e-02]),
    #     M_a = np.array([4.599896e+00,9.053326e-02,2.323585e+00]),
    #     delta_0 = np.array([3.822958e-01,5.336240e-02,6.552222e-02]),
    #     M_lna = np.array([4.423389e+00,2.202263e-01,2.165041e+00]),
    #     gamma_0 = np.array([-1.160189e+00,3.589081e-01,1.569738e-01]),
    #     M_z = np.array([-7.324986e-01,4.364703e-01,6.543076e-02]),
    #     gamma_a = np.array([-3.633671e+00,2.905324e+00,-2.849162e-01]),
    #     alpha_0 = np.array([1.973839e+00,1.305561e-01,1.766193e-02]),
    #     gamma_z = np.array([-1.218900e+00,1.134430e+00,-6.299376e-02]),
    #     alpha_a = np.array([-2.468417e+00,5.845684e-01,1.443925e+00]),
    # )
    fit_params = dict(
        epsilon_0 = np.array([-1.431495e+00,2.156496e-02,1.333165e-01]),
        alpha_lna = np.array([-1.816299e+00,4.326952e-01,1.388066e+00]),
        epsilon_a = np.array([1.757030e+00,2.236675e-01,2.916873e+00]),
        alpha_z   = np.array([1.820800e-01,2.099754e-01,7.516579e-02]),
        epsilon_lna = np.array([1.350451e+00,6.305647e-02,2.736450e+00]),
        beta_0 = np.array([4.702271e-01,6.133160e-02,3.179759e-03]),
        epsilon_z = np.array([-2.178460e-01,4.804196e-01,3.468697e-02]),
        beta_a = np.array([-8.751643e-01,7.385737e-01,2.069305e-01]),
        M_0 = np.array([1.207402e+01,1.219273e-02,8.281898e-02]),
        beta_z = np.array([-4.866420e-01,3.124524e-01,6.955197e-02]),
        M_a = np.array([4.599896e+00,9.053326e-02,2.323585e+00]),
        delta_0 = np.array([3.822958e-01,5.336240e-02,6.552222e-02]),
        M_lna = np.array([4.423389e+00,2.202263e-01,2.165041e+00]),
        gamma_0 = np.array([-1.160189e+00,3.589081e-01,1.569738e-01]),
        M_z = np.array([-7.324986e-01,4.364703e-01,6.543076e-02]),
        gamma_a = np.array([-3.633671e+00,2.905324e+00,2.849162e-01]),
        alpha_0 = np.array([1.973839e+00,1.305561e-01,1.766193e-02]),
        gamma_z = np.array([-1.218900e+00,1.134430e+00,6.299376e-02]),
        alpha_a = np.array([-2.468417e+00,5.845684e-01,1.443925e+00]),
    )
    # Compute SMHM
    M_star,delta_plus,delta_minus = np.zeros(len(mhalo)),np.zeros(len(mhalo)),np.zeros(len(mhalo))
    a = 1. / (1. + z)
    for m in range(0,len(mhalo)):
        M_star[m],delta_plus[m] = compute_smhm(fit_params,a,z,mhalo[m])
        M_star[m],delta_minus[m] = compute_smhm(fit_params,a,z,mhalo[m],err=1)

    label = r'Behroozi et al. 2019 ($z=%.2f$)'%z
    print(delta_minus,delta_plus)

    return M_star-delta_minus,M_star+delta_plus,label


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
    fig, axes = plt.subplots(1, 1,figsize=(7,7), dpi=300, facecolor='w', edgecolor='k')
    ax = axes
    ax.set_xlabel(r'$M_{\rm halo}$ [M$_{\odot}$]', fontsize=20)
    ax.set_ylabel(r'$M_{*}$ [M$_{\odot}$]', fontsize=20)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.tick_params(labelsize=16)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.minorticks_on()
    ax.tick_params(which='both',axis="both",direction="in")
    

    # Add model to plot
    minmass = np.log10(3e+9)
    maxmass = np.log10(8e11)
    ax.set_xlim([10**minmass,10**maxmass])
    ax.set_ylim([1e6,1e+11])
    mhalo = np.linspace(minmass,maxmass,100)
    mhalo = 10**mhalo
    if args.ref == 'Moster+2013':
        m_model_low,m_model_high,m_label = Moster2013(args.z[-1],mhalo)

        ax.fill_between(x=mhalo,y1=m_model_low,y2=m_model_high,color='b',alpha=0.2)
        ax.plot(mhalo[m_model_low>=1e8],m_model_low[m_model_low>=1e8],color='b')
        ax.plot(mhalo[m_model_high>=1e8],m_model_high[m_model_high>=1e8],color='b')
        ax.plot(mhalo[m_model_low<1e8],m_model_low[m_model_low<1e8],color='b',linestyle='--')
        ax.plot(mhalo[m_model_high<1e8],m_model_high[m_model_high<1e8],color='b',linestyle='--')
        ind = mhalo.searchsorted(1.5e+10)
        ax.text(1.5e+10,0.25*(m_model_high[ind]+m_model_low[ind]), m_label, fontsize=16, color='b',
                rotation=45,horizontalalignment='center', verticalalignment='center')
    if args.ref == 'Moster+2018':
        m_model_low,m_model_high,m_label = Moster2018(args.z[-1],mhalo)

        ax.fill_between(x=mhalo,y1=m_model_low,y2=m_model_high,color='b',alpha=0.2)
        ax.plot(mhalo[m_model_low>=1e8],m_model_low[m_model_low>=1e8],color='b')
        ax.plot(mhalo[m_model_high>=1e8],m_model_high[m_model_high>=1e8],color='b')
        ax.plot(mhalo[m_model_low<1e8],m_model_low[m_model_low<1e8],color='b',linestyle='--')
        ax.plot(mhalo[m_model_high<1e8],m_model_high[m_model_high<1e8],color='b',linestyle='--')
        ind = mhalo.searchsorted(1.5e+10)
        ax.text(1.5e+10,0.25*(m_model_high[ind]+m_model_low[ind]), m_label, fontsize=16, color='b',
                rotation=45,horizontalalignment='center', verticalalignment='center')
    elif args.ref == 'Behroozi+2019':
        m_model_low,m_model_high,m_label = Behroozi2019(args.z[-1],mhalo)

        ax.fill_between(x=mhalo,y1=m_model_low,y2=m_model_high,color='orange',alpha=0.2)
        ax.plot(mhalo[m_model_low>=1e8],m_model_low[m_model_low>=1e8],color='orange')
        ax.plot(mhalo[m_model_high>=1e8],m_model_high[m_model_high>=1e8],color='orange')
        ax.plot(mhalo[m_model_low<1e8],m_model_low[m_model_low<1e8],color='orange',linestyle='--')
        ax.plot(mhalo[m_model_high<1e8],m_model_high[m_model_high<1e8],color='orange',linestyle='--')
        ind = mhalo.searchsorted(8e+10)
        ax.text(8e+10,0.4*(m_model_high[ind]+m_model_low[ind]), m_label, fontsize=16, color='orange',
                rotation=40,horizontalalignment='center', verticalalignment='center')
    elif args.ref == 'Behroozi+2013':
        m_model_low,m_model_high,m_label = Behroozi2013(args.z[-1],mhalo)

        ax.fill_between(x=mhalo,y1=m_model_low,y2=m_model_high,color='orange',alpha=0.2)
        ax.plot(mhalo[m_model_low>=1e8],m_model_low[m_model_low>=1e8],color='orange')
        ax.plot(mhalo[m_model_high>=1e8],m_model_high[m_model_high>=1e8],color='orange')
        ax.plot(mhalo[m_model_low<1e8],m_model_low[m_model_low<1e8],color='orange',linestyle='--')
        ax.plot(mhalo[m_model_high<1e8],m_model_high[m_model_high<1e8],color='orange',linestyle='--')
        ind = mhalo.searchsorted(8e+10)
        ax.text(8e+10,0.4*(m_model_high[ind]+m_model_low[ind]), m_label, fontsize=16, color='orange',
                rotation=40,horizontalalignment='center', verticalalignment='center')
    elif args.ref == 'all':
        m_model_low,m_model_high,m_label = Moster2013(args.z[-1],mhalo)

        ax.fill_between(x=mhalo,y1=m_model_low,y2=m_model_high,color='b',alpha=0.2)
        ax.plot(mhalo[m_model_low>=1e8],m_model_low[m_model_low>=1e8],color='b')
        ax.plot(mhalo[m_model_high>=1e8],m_model_high[m_model_high>=1e8],color='b')
        ax.plot(mhalo[m_model_low<1e8],m_model_low[m_model_low<1e8],color='b',linestyle='--')
        ax.plot(mhalo[m_model_high<1e8],m_model_high[m_model_high<1e8],color='b',linestyle='--')
        ind = mhalo.searchsorted(6e+10)
        ax.text(7e+10,0.25*(m_model_high[ind]+m_model_low[ind]), m_label, fontsize=20, color='b',
                rotation=45,horizontalalignment='center', verticalalignment='center')
            
        m_model_low,m_model_high,m_label = Behroozi2013(args.z[-1],mhalo)

        ax.fill_between(x=mhalo,y1=m_model_low,y2=m_model_high,color='orange',alpha=0.2)
        ax.plot(mhalo[m_model_low>=1e8],m_model_low[m_model_low>=1e8],color='orange')
        ax.plot(mhalo[m_model_high>=1e8],m_model_high[m_model_high>=1e8],color='orange')
        ax.plot(mhalo[m_model_low<1e8],m_model_low[m_model_low<1e8],color='orange',linestyle='--')
        ax.plot(mhalo[m_model_high<1e8],m_model_high[m_model_high<1e8],color='orange',linestyle='--')
        ind = mhalo.searchsorted(6e+9)
        ax.text(1.2e+10,1.35*(m_model_high[ind]+m_model_low[ind]), m_label, fontsize=20, color='orange',
                rotation=39,horizontalalignment='center', verticalalignment='center')
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
            mhalo[j] = sim.galaxies[progind].halo.mass['total'].to('Msun') + sim.galaxies[progind].mass['baryon'].to('Msun') +sim.galaxies[progind].mass['halo_gas'].to('Msun')
        if not sim.simulation.physics['cr'] and not sim.simulation.physics['magnetic']:
            ax.plot(mhalo,mstellar, markersize=10, marker="d", markeredgecolor='b',
                    markerfacecolor='b',label=args.model[i],color='b')
        elif not sim.simulation.physics['cr'] and sim.simulation.physics['magnetic']:
            ax.plot(mhalo,mstellar, markersize=10,marker='*', markeredgecolor='m',
                    markerfacecolor='m',label=args.model[i],color='m')
        elif sim.simulation.physics['cr']:
            ax.plot(mhalo,mstellar, markersize=10,marker='v', markeredgecolor='g',
                    markerfacecolor='g',label=args.model[i],color='g')
        # if i == 0:
        #     # Add fixed fb to plot
        #     omega_dm = 1.0 - (sim.simulation.omega_lambda + sim.simulation.omega_k + sim.simulation.omega_baryon)
        #     fb = sim.simulation.omega_baryon / (omega_dm)
        #     print('BARYON FRACTION (Omega Baryon / Omega DM) IN SIMULATION: ',fb)
        #     mhalo_model = np.linspace(minmass,maxmass,100)
        #     mhalo_model = 10**mhalo_model
        #     mstar = lambda x: fb*x
        #     ax.plot(mhalo_model,mstar(mhalo_model),linestyle='dashdot',color='k')
        #     midx = 0.5*(maxmass+minmass)
        #     an = RotationAwareAnnotation(r'$M_*=f_{\rm b}M_{\rm DM}$', xy=(10**midx,mstar(1.15*10**midx)), p=(10**(midx+0.5),mstar(1.1*10**(midx+0.55))), ax=ax,
        #                                 xytext=(-1,1), textcoords="offset points", 
        #                                 ha="center", va="baseline", fontsize=20)
        # Compute conversion fractions just to print to screen (may be useful)
        # print(20*'-')
        # print('BARYON CONVERSION EFFICIENCY [%]')
        # print('Measured for model %s at the redshift bins of '%(args.model[i]),args.z)
        # max_mstellar = mstar(mhalo)
        # efficiency = 100 * mstellar / max_mstellar
        # print(efficiency)
        # print(20*'-')

    # Plot vertical lines separating redshifts
    ax.text(4.7e+9,2e+6, r'$z=8$', fontsize=20, color='k',alpha=0.5)
    ax.plot([1e+10,1e+10],[1e6,1e+11],'k--',alpha=0.3) # z=8
    ax.text(2e+10,2e+6, r'$z=6$', fontsize=20, color='k',alpha=0.5)
    ax.plot([5e+10,5e+10],[1e6,1e+11],'k--',alpha=0.3) # z=6
    ax.text(6e+10,2e+6, r'$z=3$', fontsize=20, color='k',alpha=0.5)
    # ax.plot([4e+9,1e+10],[1e6,1e+11],'k--',alpha=0.3) # z=3
    ax.plot([1.2e+11,1.2e+11],[1e6,1e+11],'k--',alpha=0.3) # z=1.5
    ax.text(1.7e+11,2e+6, r'$z=1.5$', fontsize=20, color='k',alpha=0.5)

    fig.subplots_adjust(top=0.97, bottom=0.1,right=0.99,hspace=0.0)
    if args.NUT:
        args.ind = 'NUT'
    fig.savefig(os.getcwd()+'/smf_'+str(args.ref)+'_'+str(args.ind)+'.pdf', format='pdf', dpi=300)
