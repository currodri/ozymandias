import numpy as np
# TODO: It should use a global fortran module, not the particles one
import sys
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/part')
import part2

class SimulationAttributes(object):
    """Class that contains the attributes of the simulation."""

    def __init__(self):
        pass
    
    def assign_attributes(self, obj):
        """Get attributes from yT and assign them to OZY object."""
        # TODO: Make sure to add an option that does not use yT but just
        # reads in the info_xxxxx.txt file.
        ds = obj.yt_dataset

        self.redshift        = ds.current_redshift
        self.scale_factor    = 1. / (1. + self.redshift)
        self.time            = ds.current_time
        self.omega_matter    = ds.omega_matter
        self.omega_lambda    = ds.omega_lambda
        self.fullpath        = ds.fullpath
        self.hubble_constant = ds.hubble_constant
        self.parameters      = ds.parameters
        self.boxsize = ds.domain_width[0].to(obj.units['length'])
        # TODO: Read this from simulation file. Right now for NUT:
        self.zoom = True
        
        H0 = ds.quan(self.hubble_constant * 100. * 3.24077929e-20, '1/s')
        Om_0 = ds.cosmology.omega_matter
        Ol_0 = ds.cosmology.omega_lambda
        Ok_0 = ds.cosmology.omega_curvature
        Or_0 = ds.cosmology.omega_radiation
        self.E_z = np.sqrt(
            Ol_0 +
            Ok_0 * (1. + self.redshift)**2. +
            Om_0 * (1. + self.redshift)**3. +
            Or_0 * (1. + self.redshift)**4.
        )
        self.Om_z = Om_0 * (1. + self.redshift)**3./(self.E_z**2.)
        self.H_z = H0 * self.E_z
        self.G = ds.quan(4.51691362044e-39, 'kpc**3/(Msun * s**2)')  # kpc^3 / (Msun s^2)
        
        self.critical_density = ds.quan(
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
        from yt import YTArray
        hdd = hd.create_group('simulation_attributes')
        units = {}
        for k,v in self.__dict__.items():
            if isinstance(v, YTArray):
                hdd.attrs.create(k, v.d)
                units[k] = v.units
            elif isinstance(v, (int, float, bool, np.number)):
                hdd.attrs.create(k, v)
            elif isinstance(v, str):
                hdd.attrs.create(k, v.encode('utf8'))
            
        uhdd = hdd.create_group('units')
        for k, v in units.items():
            uhdd.attrs.create(k, str(v).encode('utf8'))
        
        phdd = hdd.create_group('parameters')
        for k,v in self.parameters.items():
            if k != 'namelist':
                phdd.attrs.create(k, v)

        phyhdd = hdd.create_group('physics')
        for k,v in self.physics.items():
            if k != 'namelist':
                phyhdd.attrs.create(k, v)
    
    def _unpack(self, obj, hd):
        """Get, if they exist, the simulation attributes from the input HDF5 file and assign them to OZY object."""
        if 'simulation_attributes' not in hd.keys():
            print('WARNING: Simulation attributes not found in file.')
            return
        from yt.units.yt_array import YTArray
        
        hdd = hd['simulation_attributes']
        for k,v in hdd.attrs.items():
            setattr(self, k, v)
        uhdd = hdd['units']
        for k,v in uhdd.attrs.items():
            setattr(self, k, YTArray(getattr(self, k), v, registry=obj.unit_registry))
        
        phdd = hdd['parameters']
        self.parameters = {}
        for k, v in phdd.attrs.items():
            self.parameters[k] = v
