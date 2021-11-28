import numpy as np
# TODO: It should use a global fortran module, not the particles one
import sys
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/part')
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
        self.boxsize = obj.quantity(obj._info['unit_l'],'cm').to(obj.units['length'])
        # TODO: Read this from simulation file. Right now for NUT:
        self.zoom = True
        
        H0 = obj.quantity(self.hubble_constant * 100. * 3.24077929e-20, '1/s')
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
        virial_density = (177.65287921960845 * (1. + 0.4093 * (1./self.Om_z - 1.)**0.9052) - 1.) * self.Om_z
        
        self.Densities = np.array([200 * self.critical_density.to('Msun / kpc**3').d,
                                   500 * self.critical_density.to('Msun / kpc**3').d,
                                   2500 * self.critical_density.to('Msun / kpc**3').d])

        # Determine the type of physics included in this simulation
        self.physics = {'hydro':False,
                        'metals':False,
                        'magnetic':False,
                        'cr':False,
                        'rt':False,
                        'bh':False,
                        'AGN':False,
                        'dust':False}
        varIDs = part2.io_ramses.hydroID()
        part2.io_ramses.read_hydrofile_descriptor(self.fullpath,varIDs)

        if varIDs.density != 0 and varIDs.vx != 0:
            self.physics['hydro'] = True
        if varIDs.metallicity != 0:
            self.physics['metals'] = True
        if varIDs.blx != 0:
            self.physics['magnetic'] = True
        if varIDs.cr_pressure != 0:
            self.physics['cr'] = True
        if varIDs.xhii != 0:
            self.physics['rt'] = True

        
        
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
