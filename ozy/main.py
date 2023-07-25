import numpy as np
from ozy.sim_attributes import SimulationAttributes
from ozy.utils import get_mu,get_electron_mu
from unyt import UnitRegistry,unyt_array,unyt_quantity

class OZY(object):
    """Master OZY class.
    OZY objects contain all the necessary references to halos
    and galaxies in an individual simulation snapshot.

    It can be saved as a portable, standalone HDF5 file which allows
    general analysis without requiring the original snapshot.
    """
    def __init__(self, fullpath, *args, **kwargs):
        self._args   = args
        self._kwargs = kwargs

        self.units = dict(
            mass        = 'Msun',
            length      = 'kpc',
            velocity    = 'km/s',
            time        = 'yr',
            temperature = 'K'
        )

        self._info = self._get_my_info(fullpath)
        self.unit_registry = self._get_unit_registry()
        self.simulation  = SimulationAttributes()
        self._assign_simulation_attributes(fullpath)
        
        self.nhalos      = 0
        self.ngalaxies   = 0
        self.halos       = []
        self.galaxies    = []
        self.group_types = []
        
        if 'catalogue_reader' in self._kwargs:
            self.catalogue_reader = self._kwargs['catalogue_reader']
        else:
            print('Reverting to the HaloMaker reader, since no catalogue_reader was provided!')
            from ozy.read_HaloMaker import build_HaloMaker
            self.catalogue_reader = build_HaloMaker

    def array(self, value, units):
        return unyt_array(value, units, registry=self.unit_registry)

    def quantity(self, value, units):
        return unyt_quantity(value, units, registry=self.unit_registry)

    @property
    def _has_halos(self):
        """Check if the dataset has halos."""
        if self.nhalos > 0:
            return True
        else:
            return False

    @property
    def _has_galaxies(self):
        """Check if the dataset has galaxies."""
        if self.ngalaxies > 0:
            return True
        else:
            return False

    def _get_my_info(self,fullpath):
        from ozy.utils import read_infofile
        index = int(fullpath[-5:])
        infofile_path = fullpath + '/info_%05d.txt' % index
        return read_infofile(infofile_path)

    def _get_unit_registry(self):
        from unyt.dimensions import length,mass,time,temperature,dimensionless,magnetic_field_cgs
        from unyt import mp,kb,erg,K,g
        
        registry = UnitRegistry(unit_system='cgs')

        _X = 0.76  # H fraction, hardcoded
        _Y = 0.24  # He fraction, hardcoded
        mean_molecular_weight_factor = get_mu(_X, _Y)
        electron_molecular_weight_factor = get_electron_mu(_X, _Y)

        # unyt stores internally in MKS units (m,kg,s), so a couple of
        # transformations are required to the CGS units in RAMSES
        length_unit = self._info["unit_l"]/1e+2 # cm to m
        density_unit = self._info["unit_d"] * 1e+3 # g/cm**3 to kg/m**3
        time_unit = self._info["unit_t"]
        mass_unit = density_unit * length_unit ** 3
        magnetic_unit = np.sqrt(4. *np.pi) * length_unit * (density_unit**0.5) / time_unit
        velocity_unit = length_unit / time_unit
        pressure_unit = density_unit * (length_unit / time_unit) ** 2
        temperature_unit = velocity_unit ** 2 * mp.to('kg').d * mean_molecular_weight_factor / kb.to('kg*m**2/(K*s**2)').d
        s_entropy_unit = 1.4e+8 * erg / K / g
        s_entropy_unit = float(s_entropy_unit.to('m**2/s**2/K').d)
        pseudo_entropy_unit = (mp.to('kg').d**(5./3.)) * mean_molecular_weight_factor * electron_molecular_weight_factor**(2./3.) * pressure_unit / (density_unit**(5./3.))


        # Code length
        registry.add("code_length", base_value=length_unit, dimensions=length)
        # Code time
        registry.add("code_time", base_value=time_unit, dimensions=time)
        # Code density
        registry.add("code_density", base_value=density_unit, dimensions=mass/(length**3))
        # Code mass
        registry.add("code_mass", base_value=density_unit*(length_unit**3), dimensions=mass)
        # Code velocity
        registry.add("code_velocity", base_value=velocity_unit, dimensions=length/time)
        # Code pressure
        registry.add("code_pressure", base_value=pressure_unit, dimensions=(mass)/((length)*(time)**2))
        # Code energy
        registry.add("code_energy", base_value=mass_unit*velocity_unit**2, 
                    dimensions=mass*(length**2)/(time**2))
        # Code specific energy
        registry.add("code_specific_energy", base_value=velocity_unit**2, 
                    dimensions=(length**2)/(time**2))

        # Code magnetic in Lorentz-Heavyside rational units
        registry.add("code_magnetic", base_value=magnetic_unit, 
                    dimensions=magnetic_field_cgs)
        # Code magnetic in standard cgs units
        registry.add("code_magnetic_standard", base_value=np.sqrt(4*np.pi)*magnetic_unit, 
                    dimensions=magnetic_field_cgs)
        # Code energy density
        registry.add("code_energy_density", base_value=magnetic_unit**2, 
                    dimensions=mass/(length*time**2))
        # Code temperature
        registry.add("code_temperature", base_value=temperature_unit, dimensions=temperature)
        # Code metallicity
        registry.add("code_metallicity", base_value=1.0, dimensions=dimensionless)
        # Code specific entropy
        registry.add("code_specific_entropy", base_value=s_entropy_unit, 
                    dimensions=(length**2)/(temperature*time**2))
        # Code pseudo-entropy
        registry.add('code_pseudo_entropy',base_value=pseudo_entropy_unit,
                    dimensions=mass*(length**4)/(time**2))

        return registry


    def _assign_simulation_attributes(self,fullpath):
        """Assign simulation attributes to the OZY object, if it has not been done before."""
        self.simulation.assign_attributes(self,fullpath)
    
    def _assign_groups(self):
        """Assign galaxies to halos to galaxies.
            Also connect halos with their central galaxy."""
        import ozy.group_assignment as assign
        assign.galaxies_to_halos(self)
        assign.central_galaxies(self)
    
    def _link_groups(self):
        """Two-way linking of objects."""
        from ozy.group_linking import link
        link.galaxies_to_halos(self)
        link.create_sublists(self)
    
    def save(self, filename):
        """Save OZY object as HDF5 file."""
        from ozy.saver import save
        save(self, filename)
    
    def run_halocatalogue(self, *args, **kwargs):
        """This is the central function of the OZY class for the HALOMAKER catalogues.

        This method is reponsible for:
        1) Calling the halo catalogue reader to extract the basic information
        2) Creating halos and galaxies
        3) Linking objects through the chosen method
        4) Computing additional quantities

        """
        import ozy.group_assignment as assign
        import ozy.group_linking as link
        
        self._args = args
        self._kwargs = kwargs

        # Run the provided halo catalogue reader
        self.catalogue_reader(self,*args,**kwargs)

        if self._has_halos:
            # Make assignment
            assign.galaxies_to_halos(self)

            # Now process galaxies using their assigned halo
            main_gal_ID = -1
            if 'main_gal' in self._kwargs:
                if self._kwargs['main_gal']:
                    print('Computing details just for main galaxy in simulation.')
                    masses = [i.virial_quantities['mass'] for i in self.galaxies]
                    main_gal_ID = self.galaxies[np.argmax(masses)].ID
                    
            if main_gal_ID != -1:
                for gal in self.galaxies:
                    gal._empty_galaxy()
                self.galaxies[np.argmax(masses)]._process_galaxy()
            else:
                for gal in self.galaxies:
                    gal._process_galaxy(**self._kwargs)

            # Link objects between each other
            link.galaxies_to_halos(self)

            assign.central_galaxies(self)
            link.create_sublists(self)
        else:
            print("WARNING: Not a single virialised halo above the minimum particle threshold.")

    def galaxies_summary(self, top=10):
        """Method to briefly print information for the most massive galaxies in the catalogue."""
        from ozy.utils import info_printer
        info_printer(self, 'galaxy', top)

    def halos_summary(self, top=10):
        """Method to briefly print information for the most massive halos in the catalogue."""
        from ozy.utils import info_printer
        info_printer(self, 'halo', top)

