import numpy as np


class RAMSESCooling:
    """Vectorized Python translation of the RAMSES Fortran cooling module.

    Computes H+He primordial cooling/heating and metal-line cooling rates
    identical to those used by RAMSES when haardt_madau=.true. (Courty UV
    background) and z_reion=<z_reion> are set in the namelist.

    Usage
    -----
    engine = RAMSESCooling(z_reion=9.0)
    L_rad, H_uvb = engine.get_cooling_heating(T, n_H, z, Z_solar, boost)

    Parameters
    ----------
    T         : ndarray, physical temperature [K]
    n_H       : ndarray, hydrogen number density [cm^-3]
    z         : float,   redshift
    Z_solar   : float,   metallicity in solar units
    boost     : ndarray, UV self-shielding boost factor (1 = no shielding)

    Returns
    -------
    L_rad  : volumetric cooling rate [erg cm^-3 s^-1]
    H_uvb  : volumetric UV photo-heating rate [erg cm^-3 s^-1]

    Notes
    -----
    Metal cooling replicates RAMSES cmp_metals exactly:
      metal = (Cloudy07_table + fine_structure) * f_courty(z, nH, T)
    Compton cooling/heating is not included (it lives in the separate
    cool_com/heat_com columns of the exported RAMSES table).
    """

    def __init__(self, z_reion=9.0, haardt_madau=True):
        if not haardt_madau:
            raise NotImplementedError(
                "Only the Courty/Haardt-Madau UV background "
                "(haardt_madau=True) is implemented.")
        self.haardt_madau = haardt_madau
        self.eV2erg  = 1.60217653e-12
        self.kB      = 1.380649e-16
        self.mH      = 1.6735e-24
        self.YHelium = 0.24
        self.zreioniz = z_reion

        # Courty UV background polynomial coefficients (order 0–7)
        # Rows: HI ion, HeI ion, HeII ion, HI heat, HeI heat, HeII heat
        self.coefcourty = np.array([
            [-13.5857,  1.24475,    0.187739, -0.430409,  0.152544, -0.0246448,  0.00192622, -5.89772e-05],
            [-14.0242,  1.99211,   -0.490766, -0.122646,  0.0776501,-0.0146310,  0.00123335, -3.96066e-05],
            [-15.6627,  0.128240,   1.65633,  -1.23799,   0.372157, -0.0561687,  0.00422696, -0.000126344],
            [-24.8422,  1.50750,   -0.0699428,-0.308682,  0.122196, -0.0205179,  0.00163695, -5.08050e-05],
            [-25.0252,  1.79577,   -0.159054, -0.300924,  0.125343, -0.0214598,  0.00173377, -5.43576e-05],
            [-26.4168,  0.0479454,  1.70948,  -1.26395,   0.378922, -0.0570957,  0.00428897, -0.000127909]
        ])
        self.coef_fit = np.array([20., 20., 20., 20., 20., 20.])
        self.beta_fit = np.array([6,   6,   8,   6,   6,   8])

        # Cloudy 07 metal cooling table  (log10 T [K], log10 Lambda/nH^2 [erg cm^3 s^-1])
        self.temp_cc07 = np.array([
            3.9684,4.0187,4.0690,4.1194,4.1697,4.2200,4.2703,4.3206,4.3709,4.4212,
            4.4716,4.5219,4.5722,4.6225,4.6728,4.7231,4.7734,4.8238,4.8741,4.9244,
            4.9747,5.0250,5.0753,5.1256,5.1760,5.2263,5.2766,5.3269,5.3772,5.4275,
            5.4778,5.5282,5.5785,5.6288,5.6791,5.7294,5.7797,5.8300,5.8804,5.9307,
            5.9810,6.0313,6.0816,6.1319,6.1822,6.2326,6.2829,6.3332,6.3835,6.4338,
            6.4841,6.5345,6.5848,6.6351,6.6854,6.7357,6.7860,6.8363,6.8867,6.9370,
            6.9873,7.0376,7.0879,7.1382,7.1885,7.2388,7.2892,7.3395,7.3898,7.4401,
            7.4904,7.5407,7.5911,7.6414,7.6917,7.7420,7.7923,7.8426,7.8929,7.9433,
            7.9936,8.0439,8.0942,8.1445,8.1948,8.2451,8.2955,8.3458,8.3961,8.4464,8.4967
        ])
        self.excess_cooling_cc07 = np.array([
            -24.9082,-24.9082,-24.5503,-24.0898,-23.5328,-23.0696,-22.7758,-22.6175,
            -22.5266,-22.4379,-22.3371,-22.2289,-22.1181,-22.0078,-21.8992,-21.7937,
            -21.6921,-21.5961,-21.5089,-21.4343,-21.3765,-21.3431,-21.3274,-21.3205,
            -21.3142,-21.3040,-21.2900,-21.2773,-21.2791,-21.3181,-21.4006,-21.5045,
            -21.6059,-21.6676,-21.6877,-21.6934,-21.7089,-21.7307,-21.7511,-21.7618,
            -21.7572,-21.7532,-21.7668,-21.7860,-21.8129,-21.8497,-21.9035,-21.9697,
            -22.0497,-22.1327,-22.2220,-22.3057,-22.3850,-22.4467,-22.4939,-22.5205,
            -22.5358,-22.5391,-22.5408,-22.5408,-22.5475,-22.5589,-22.5813,-22.6122,
            -22.6576,-22.7137,-22.7838,-22.8583,-22.9348,-23.0006,-23.0547,-23.0886,
            -23.1101,-23.1139,-23.1147,-23.1048,-23.1017,-23.0928,-23.0969,-23.0968,
            -23.1105,-23.1191,-23.1388,-23.1517,-23.1717,-23.1837,-23.1986,-23.2058,
            -23.2134,-23.2139,-23.2107
        ])

        # Courty UV flux table φ(z) used by f_courty
        self.z_courty = np.array([
            0.00000, 0.04912, 0.10060, 0.15470, 0.21140, 0.27090, 0.33330, 0.39880,
            0.46750, 0.53960, 0.61520, 0.69450, 0.77780, 0.86510, 0.95670, 1.05300,
            1.15400, 1.25900, 1.37000, 1.48700, 1.60900, 1.73700, 1.87100, 2.01300,
            2.16000, 2.31600, 2.47900, 2.64900, 2.82900, 3.01700, 3.21400, 3.42100,
            3.63800, 3.86600, 4.10500, 4.35600, 4.61900, 4.89500, 5.18400, 5.48800,
            5.80700, 6.14100, 6.49200, 6.85900, 7.24600, 7.65000, 8.07500, 8.52100,
            8.98900, 9.50000
        ])
        self.phi_courty = np.array([
            0.0499886, 0.0582622, 0.0678333, 0.0788739, 0.0915889, 0.1061913, 0.1229119,
            0.1419961, 0.1637082, 0.1883230, 0.2161014, 0.2473183, 0.2822266, 0.3210551,
            0.3639784, 0.4111301, 0.4623273, 0.5172858, 0.5752659, 0.6351540, 0.6950232,
            0.7529284, 0.8063160, 0.8520859, 0.8920522, 0.9305764, 0.9682031, 1.0058810,
            1.0444020, 1.0848160, 1.1282190, 1.1745120, 1.2226670, 1.2723200, 1.3231350,
            1.3743020, 1.4247480, 1.4730590, 1.5174060, 1.5552610, 1.5833640, 1.5976390,
            1.5925270, 1.5613110, 1.4949610, 1.3813710, 1.2041510, 0.9403100, 0.5555344,
            0.0000000
        ])

    # ------------------------------------------------------------------
    # UV self-shielding suppression of metal cooling
    # ------------------------------------------------------------------

    def get_f_courty(self, z, nH, T):
        """UV self-shielding suppression factor for metal cooling.

        Implements RAMSES cmp_metals: f = 1/(1 + ux/g)
          ux = 1e-4 * phi(z) / nH
          g  = 0.4*(T/1e5)^0.15 + 10*exp(-1e6/T)

        Uses the exact Fortran index formula (uniform-grid approximation on the
        non-uniform z_courty array) to match the RAMSES table bit-for-bit.
        """
        c1, c2, TT0, TTC, alpha1 = 0.4, 10.0, 1e5, 1e6, 0.15
        if z <= 0.0 or z >= self.z_courty[-1]:
            phi_z = 0.0
        else:
            iz = int(z / self.z_courty[-1] * 49.)
            iz = min(iz, 48)
            iz = max(iz, 0)
            dz = self.z_courty[iz + 1] - self.z_courty[iz]
            phi_z = (self.phi_courty[iz + 1] * (z - self.z_courty[iz])
                     + self.phi_courty[iz]   * (self.z_courty[iz + 1] - z)) / dz
        ux = 1e-4 * phi_z / np.maximum(nH, 1e-300)
        g  = c1 * (T / TT0)**alpha1 + c2 * np.exp(-TTC / np.maximum(T, 1.0))
        return 1.0 / (1.0 + ux / np.maximum(g, 1e-300))

    # ------------------------------------------------------------------
    # UV photo-ionization and photo-heating rates (Courty model)
    # ------------------------------------------------------------------

    def get_courty_rates(self, z):
        """Courty UV background rates at redshift z.

        Returns array of 6 rates: [HI_ion, HeI_ion, HeII_ion,
                                    HI_heat, HeI_heat, HeII_heat]
        All rates are zero for z >= z_reion.
        """
        if z >= self.zreioniz:
            return np.zeros(6)
        zz = max(z, 1.0e-15)
        rates = np.zeros(6)
        for i in range(6):
            hh = sum(self.coefcourty[i, j] * (zz**j) for j in range(8))
            hhreion = self.coef_fit[i] * (zz / self.zreioniz)**self.beta_fit[i]
            rates[i] = 10.0**(hh - hhreion)
        return rates

    # ------------------------------------------------------------------
    # Primordial reaction rate coefficients (vectorized over T arrays)
    # ------------------------------------------------------------------

    def cool_bre(self, ispec, T):
        term = 1.1 + 0.34 * np.exp(-(5.5 - np.log10(T))**2 / 3.0)
        if ispec == 0: return 1.42e-27 * np.sqrt(T) * term
        if ispec == 1: return 1.42e-27 * np.sqrt(T) * term
        if ispec == 2: return 5.68e-27 * np.sqrt(T) * term

    def cool_exc(self, ispec, T):
        T5 = 1e-5 * T
        if ispec == 0: return 7.50e-19 / (1. + np.sqrt(T5)) * np.exp(-118348. / T)
        if ispec == 1: return 9.10e-27 / (1. + np.sqrt(T5)) / (T**0.1687) * np.exp(-13179. / T)
        if ispec == 2: return 5.54e-17 / (1. + np.sqrt(T5)) / (T**0.397)  * np.exp(-473638. / T)

    def cool_rec(self, ispec, T):
        T3 = 1e-3 * T; T6 = 1e-6 * T
        if ispec == 0: return 8.70e-27 * np.sqrt(T) / (T3**0.2) / (1. + T6**0.7)
        if ispec == 1: return 1.55e-26 * T**0.3647
        if ispec == 2: return 3.48e-26 * np.sqrt(T) / (T3**0.2) / (1. + T6**0.7)

    def taux_die(self, T):
        return 1.9e-3 * T**(-1.5) * np.exp(-470000. / T) * (1. + 0.3 * np.exp(-94000. / T))

    def cool_die(self, T):
        return 1.24e-13 * T**(-1.5) * np.exp(-470000. / T) * (1. + 0.3 * np.exp(-94000. / T))

    def taux_rec(self, ispec, T):
        T3 = 1e-3 * T; T6 = 1e-6 * T
        if ispec == 0: return 0.75 * 8.40e-11 / np.sqrt(T) / (T3**0.2) / (1. + T6**0.7)
        if ispec == 1: return 1.50e-10 / T**0.6353 + self.taux_die(T)
        if ispec == 2: return 3.36e-10 / np.sqrt(T) / (T3**0.2) / (1. + T6**0.7)

    def taux_ion(self, ispec, T):
        T5 = 1e-5 * T
        if ispec == 0: return 2.0 * 5.85e-11 * np.sqrt(T) / (1. + np.sqrt(T5)) * np.exp(-157809.1 / T)
        if ispec == 1: return 2.0 * 2.38e-11 * np.sqrt(T) / (1. + np.sqrt(T5)) * np.exp(-285335.4 / T)
        if ispec == 2: return 2.0 * 5.68e-12 * np.sqrt(T) / (1. + np.sqrt(T5)) * np.exp(-631515.0 / T)

    def cool_ion(self, ispec, T):
        T5 = 1e-5 * T
        if ispec == 0: return 2.0 * 1.27e-21 * np.sqrt(T) / (1. + np.sqrt(T5)) * np.exp(-157809.1 / T)
        if ispec == 1: return 2.0 * 9.38e-22 * np.sqrt(T) / (1. + np.sqrt(T5)) * np.exp(-285335.4 / T)
        if ispec == 2: return 2.0 * 4.95e-22 * np.sqrt(T) / (1. + np.sqrt(T5)) * np.exp(-631515.0 / T)

    # ------------------------------------------------------------------
    # Chemical equilibrium solver
    # ------------------------------------------------------------------

    def cmp_chem_eq(self, T, n_H, t_rad_spec):
        """Iterative H/He ionization equilibrium (25 iterations).

        Returns (n_E, n_HI, n_HII, n_HEI, n_HEII, n_HEIII).
        """
        yy = self.YHelium / (1. - self.YHelium) / 4.
        t_rad_HI, t_rad_HEI, t_rad_HEII = t_rad_spec
        t_rec_HI   = self.taux_rec(0, T);  t_ion_HI   = self.taux_ion(0, T)
        t_rec_HEI  = self.taux_rec(1, T);  t_ion_HEI  = self.taux_ion(1, T)
        t_rec_HEII = self.taux_rec(2, T);  t_ion_HEII = self.taux_ion(2, T)
        n_E = np.copy(n_H)
        for _ in range(25):
            n_E_safe    = np.maximum(n_E, 1e-15 * n_H)
            t_ion2_HI   = t_ion_HI   + t_rad_HI   / n_E_safe
            t_ion2_HEI  = t_ion_HEI  + t_rad_HEI  / n_E_safe
            t_ion2_HEII = t_ion_HEII + t_rad_HEII / n_E_safe
            n_HI  = t_rec_HI  / (t_ion2_HI  + t_rec_HI)  * n_H
            n_HII = t_ion2_HI / (t_ion2_HI  + t_rec_HI)  * n_H
            x1 = (t_rec_HEII * t_rec_HEI
                  + t_ion2_HEI  * t_rec_HEII
                  + t_ion2_HEII * t_ion2_HEI)
            n_HEIII = yy * t_ion2_HEII * t_ion2_HEI / x1 * n_H
            n_HEII  = yy * t_ion2_HEI  * t_rec_HEII / x1 * n_H
            n_HEI   = yy * t_rec_HEII  * t_rec_HEI  / x1 * n_H
            n_E = 0.5 * n_E + 0.5 * (n_HII + n_HEII + 2.0 * n_HEIII)
        return n_E, n_HI, n_HII, n_HEI, n_HEII, n_HEIII

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def get_cooling_heating(self, T, n_H, z, Z_solar, boost_array):
        """Total volumetric cooling and UV heating rates.

        Parameters
        ----------
        T          : array, physical temperature [K]
        n_H        : array, H number density [cm^-3]
        z          : float, redshift
        Z_solar    : float, metallicity in solar units
        boost_array: array, UV self-shielding attenuation (1 = none)

        Returns
        -------
        L_rad : array, [erg cm^-3 s^-1]  — radiative cooling
        H_uvb : array, [erg cm^-3 s^-1]  — UV photo-heating
        """
        uv = self.get_courty_rates(z)
        t_rad_spec = np.array([uv[0], uv[1], uv[2]])[:, None] * boost_array
        h_rad_spec = np.array([uv[3], uv[4], uv[5]])[:, None] * boost_array

        n_E, n_HI, n_HII, n_HEI, n_HEII, n_HEIII = self.cmp_chem_eq(T, n_H, t_rad_spec)

        prim = (
            self.cool_bre(0, T) * n_E * n_HII   / n_H**2
          + self.cool_bre(1, T) * n_E * n_HEII  / n_H**2
          + self.cool_bre(2, T) * n_E * n_HEIII / n_H**2
          + self.cool_ion(0, T) * n_E * n_HI    / n_H**2
          + self.cool_ion(1, T) * n_E * n_HEI   / n_H**2
          + self.cool_ion(2, T) * n_E * n_HEII  / n_H**2
          + self.cool_rec(0, T) * n_E * n_HII   / n_H**2
          + self.cool_rec(1, T) * n_E * n_HEII  / n_H**2
          + self.cool_rec(2, T) * n_E * n_HEIII / n_H**2
          + self.cool_die(T)    * n_E * n_HEII  / n_H**2
          + self.cool_exc(0, T) * n_E * n_HI    / n_H**2
          + self.cool_exc(1, T) * n_E * n_HEI   / n_H**2
          + self.cool_exc(2, T) * n_E * n_HEII  / n_H**2
        )

        # Metal cooling: Cloudy 07 + fine-structure infrared, both × f_courty
        # Replicates RAMSES cmp_metals: metal_tot = (10^lcool1 + 10^lcool2) * f_courty
        T_guard    = np.maximum(T, 10.)
        T_max_cc07 = 10.**self.temp_cc07[-1]      # ≈ 3.14e8 K
        logT       = np.log10(T_guard)
        logT_c     = np.clip(logT, self.temp_cc07[0], self.temp_cc07[-1])
        metal_cc07 = np.where(
            (T_guard >= 10.**self.temp_cc07[0]) & (T_guard < T_max_cc07),
            10.**np.interp(logT_c, self.temp_cc07, self.excess_cooling_cc07), 0.0)
        metal_fs = np.where(
            T_guard < T_max_cc07,
            10.**(-31.522879 + 2.0*logT - 20.0/T_guard - T_guard*4.342944e-5), 0.0)
        f_cty = self.get_f_courty(z, n_H, T_guard)
        metal = (metal_cc07 + metal_fs) * f_cty * Z_solar

        uvb_heat = (h_rad_spec[0] * n_HI + h_rad_spec[1] * n_HEI
                    + h_rad_spec[2] * n_HEII) / n_H**2

        return (prim + metal) * n_H**2, uvb_heat * n_H**2
