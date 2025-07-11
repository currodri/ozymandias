import numpy as np
import part2

class SimulationAttributes(object):
    """Class that contains the attributes of the simulation."""

    def __init__(self):
        pass
   
    def assign_attributes(self, obj,fullpath):
        """Get attributes from info file and assign them to OZY object."""

        self.redshift        = obj._info['redshift']
        self.scale_factor    = obj._info['aexp']
        self.time            = obj._info['time']
        self.omega_matter    = obj._info['omega_m']
        self.omega_lambda    = obj._info['omega_l']
        self.omega_k         = obj._info['omega_k']
        self.omega_baryon    = obj._info['omega_b']
        self.fullpath        = fullpath
        self.hubble_constant = obj._info['H0']
        self.boxsize = obj.quantity(obj._info['unit_l']*obj._info['boxlen'],'cm').to(obj.units['length'])
        # Get from the DM particle numbers whether or not this is a zoom simulation
        self.zoom = False
        if self.omega_matter!=1.0:
            try:
                ndm = obj._header['particle_counts']['dark_matter_particles']
            except:
                ndm = obj._header['particle_counts']['DM']
            if ndm != (2**obj._info['levelmin'])**obj._info['ndim']:
                self.zoom = True
        
        H0 = obj.quantity(self.hubble_constant * 3.24077929e-20, '1/s')
        Om_0 = obj._info['omega_m']
        Ol_0 = obj._info['omega_l']
        Ok_0 = obj._info['omega_k']
        Or_0 = 0.0 # TODO Figure out why there's not omega_r in RAMSES info file
        self.E_z = np.sqrt(
            Ol_0 +
            Ok_0 * (1. + self.redshift)**2. +
            Om_0 * (1. + self.redshift)**3. +
            Or_0 * (1. + self.redshift)**4.
        )
        self.Om_z = Om_0 * (1. + self.redshift)**3./(self.E_z**2.)
        self.H_z = H0 * self.E_z
        self.G = obj.quantity(4.51691362044e-39, 'kpc**3/(Msun * s**2)')  # kpc^3 / (Msun s^2)
        
        self.critical_density = obj.quantity(
            (3. * self.H_z**2) / (8. * np.pi * self.G.d),
            'Msun / kpc**3'
        )
        self.virial_density = (177.65287921960845 * (1. + 0.4093 * (1./self.Om_z - 1.)**0.9052) - 1.) * self.Om_z
        
        self.Densities = np.array([200 * self.critical_density.to('Msun / kpc**3').d,
                                   500 * self.critical_density.to('Msun / kpc**3').d,
                                   2500 * self.critical_density.to('Msun / kpc**3').d])

        # Determine the type of physics included in this simulation
        self.physics = {'hydro':False,
                        'metallicity':False,
                        'metals':False,
                        'CO':False,
                        'magnetic':False,
                        'cr':False,
                        'rt':False,
                        'bh':False,
                        'AGN':False,
                        'dust':False,
                        'PAHs':False}
        
        if obj.vardict.get('density') !=0 and obj.vardict.get('velocity_x') != 0:
            self.physics['hydro'] = True
        if obj.vardict.get('metallicity') != 0:
            self.physics['metallicity'] = True
        if obj.vardict.get('iron_fraction') != 0:
            self.physics['metals'] = True
        if obj.vardict.get('B_left_x') != 0:
            self.physics['magnetic'] = True
        if obj.vardict.get('cr_pressure') != 0:
            self.physics['cr'] = True
        if obj.vardict.get('xHII') != 0:
            self.physics['rt'] = True
        if obj.vardict.get('CSmall') != 0:
            self.physics['dust'] = True
        if obj.vardict.get('PAHSmall') or obj.vardict.get('PAH') != 0:
            self.physics['PAHs'] = True
        if obj.vardict.get('CO_fraction') != 0:
            self.physics['CO'] = True

        # Determine the type of numerics included in this simulation
        self.numerics = {'families':False}

        if obj.part_vardict.get('family') != 0:
            self.numerics['families'] = True

        # Determine what elements are included in this simulation
        if self.physics['metals']:
            self.elements = {'helium':False,
                             'carbon':False,
                             'nitrogen':False,
                             'oxygen':False,
                             'neon':False,
                             'magnesium':False,
                             'silicon':False,
                             'sulfur':False,
                             'iron':False,
                             'calcium':False,
                             'fluoride':False}
            if obj.vardict.get('helium_fraction') != 0:
                self.elements['helium'] = True
            if obj.vardict.get('carbon_fraction') != 0:
                self.elements['carbon'] = True
            if obj.vardict.get('nitrogen_fraction') != 0:
                self.elements['nitrogen'] = True
            if obj.vardict.get('oxygen_fraction') != 0:
                self.elements['oxygen'] = True
            if obj.vardict.get('neon_fraction') != 0:
                self.elements['neon'] = True
            if obj.vardict.get('magnesium_fraction') != 0:
                self.elements['magnesium'] = True
            if obj.vardict.get('silicon_fraction') != 0:
                self.elements['silicon'] = True
            if obj.vardict.get('sulfur_fraction') != 0:
                self.elements['sulfur'] = True
            if obj.vardict.get('iron_fraction') != 0:
                self.elements['iron'] = True
            if obj.vardict.get('calcium_fraction') != 0:
                self.elements['calcium'] = True
            if obj.vardict.get('fluoride_fraction') != 0:
                self.elements['fluoride'] = True
        
        
    def _serialise(self, obj, hd):
        """This makes possible to save the simulation attributes as dataset attributes of an HDF5 file."""
        from unyt import unyt_quantity
        hdd = hd.create_group('simulation_attributes')
        units = {}
        for k,v in self.__dict__.items():
            if isinstance(v, unyt_quantity):
                hdd.attrs.create(k, v.d)
                units[k] = v.units
            elif isinstance(v, (int, float, bool, np.number)):
                hdd.attrs.create(k, v)
            elif isinstance(v, str):
                hdd.attrs.create(k, v.encode('utf8'))
            
        uhdd = hdd.create_group('units')
        for k, v in units.items():
            uhdd.attrs.create(k, str(v).encode('utf8'))

        phyhdd = hdd.create_group('physics')
        for k,v in self.physics.items():
            if k != 'namelist':
                phyhdd.attrs.create(k, v)
    
    def _unpack(self, obj, hd):
        """Get, if they exist, the simulation attributes from the input HDF5 file and assign them to OZY object."""
        if 'simulation_attributes' not in hd.keys():
            print('WARNING: Simulation attributes not found in file.')
            return

        hdd = hd['simulation_attributes']
        for k,v in hdd.attrs.items():
            setattr(self, k, v)
        uhdd = hdd['units']
        for k,v in uhdd.attrs.items():
            setattr(self, k, obj.quantity(getattr(self, k), v))
        
        try:
            phyhdd = hdd['physics']
            self.physics = {'hydro':False,
                            'metals':False,
                            'magnetic':False,
                            'cr':False,
                            'rt':False,
                            'bh':False,
                            'AGN':False,
                            'dust':False}
            for k,v in phyhdd.attrs.items():
                self.physics[k] = v
        except:
            print('No physics details in this OZY file!')
