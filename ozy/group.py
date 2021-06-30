import numpy as np
from pprint import pprint
import sys
sys.path.append('/mnt/extraspace/currodri/Codes/ozymandias/ozy/amr')
sys.path.append('/mnt/extraspace/currodri/Codes/ozymandias/ozy/part')
from amr2 import vectors
from amr2 import geometrical_regions as geo
from amr2 import filtering
from amr2 import amr_integrator
from part2 import part_integrator

from ozy.profiles import init_region

MINIMUM_STARS_PER_GALAXY = 100
MINIMUM_DM_PER_HALO      = 0

grouptypes = dict(
    halo='halos',
    galaxy='galaxies',
)

info_blacklist = [
    'obj', 'halo', 'galaxies', 'satellites',
    'galaxy_index_list_end', 'galaxy_index_list_start']

# class GroupList(object):
#     """Class to hold particle index lists.
#     """
#     def __init__(self, name):
#         self.name = name

#     def __get__(self, instance, owner):
#         if not hasattr(instance, '_%s' % self.name) or isinstance(getattr(instance, '_%s' % self.name), int):
#             from ozy.loader import restore_single_list
#             restore_single_list(instance.obj, instance, self.name)
#         return getattr(instance, instance, self.name)
    
#     def __set__(self, instance, value):
#         setattr(instance, '_%s' % self.name, value)

class Group(object):
    """This is the parent class from which the rest of groups are derived."""

    # slist = GroupList('slist')
    
    def __init__(self, obj):
        """Initialise basic properties of group object."""
        self.obj      = obj
        self.ID       = -1
        self.npart    = 0
        self.units    = 'code'
        self.position = np.array([-1, -1, -1])
        self.radius   = {}
        self.shape = {}
        self.mass     = {}
        self.velocity = np.array([0, 0, 0])
        self.angular_mom = {}
        self.virial_quantities = {}
        self.energies = {}
    
    @property
    def _valid(self):
        """Check against the minimum number of particles to see if
        this object is 'valid'."""
        if self.obj_type == 'halo' and self.ndm < MINIMUM_DM_PER_HALO:
            return False
        elif self.obj_type == 'galaxy' and self.nstar < MINIMUM_STARS_PER_GALAXY:
            return False
        else:
            return True

    def _info(self):
        """Method to quickly print out object attributes."""
        attrdict = {}
        for k,v in six.iteritems(self.__dict__):
            if k in info_blacklist:
                continue
            attrdict[k] = v
        pprint(attrdict)
        attrdict = None

class Galaxy(Group):
    """Galaxy class which has the central boolean."""
    obj_type = 'galaxy'
    def __init__(self, obj):
        super(Galaxy, self).__init__(obj)
        self.central = False
        self.halo    = None
        self.gas_density = {}
        self.sfr     = {}
        self.energies = {}
        self.metallicity = {}
        self.temperature = {}
        self.radiation = {}
    def _process_galaxy(self):
        """Process each galaxy after creation. This means
        calculating the total mass, and then calculate the rest of masses,
        radial quantities, velocity dispersions, angular momentum...
        """
        self._calculate_stardm_quantities()
        self._calculate_velocity_dispersions()
        self._calculate_gas_quantities()

        self.mass['baryon'] = self.mass['stellar'] + self.mass['gas']

        self.gas_fraction = 0.0
        if self.mass['baryon'] > 0:
            self.gas_fraction = self.mass['gas'].d / self.mass['baryon'].d

    def _calculate_stardm_quantities(self):
        """Calculate total stellar and DM masses as well as
            some useful global quantities.
        """
        output_path = self.obj.simulation.fullpath

        # Initialise region
        selected_reg = init_region(self,'sphere')

        # We do not want any particular filter, just simple integration will do
        filt = filtering.filter()

        # Initialise Fortran derived type with attributes
        # This object hold the following attributes:
        # - nstar: number of star particles in region
        # - ndm:  number of DM particles in region
        # - nvars: number of variables
        # - nwvars:  number of variables for weighting
        # - varnames: names of variables
        # - wvarnames:  names of variables for weighting
        # - data: organised in numpy array of shape (nvars,nwvars,4)
        #           each of those 4 values are (final, min, max, sum of weights)
        glob_attrs = part_integrator.part_region_attrs()
        glob_attrs.nvars = 5
        glob_attrs.nwvars = 2
        part_integrator.allocate_part_regions_attrs(glob_attrs)
        glob_attrs.varnames.T.view('S128')[0] = b'star/mass'.ljust(128)
        # TODO: SFR indicators should be taken as arguments
        glob_attrs.varnames.T.view('S128')[1] = b'star/sfr_10'.ljust(128)
        glob_attrs.varnames.T.view('S128')[2] = b'star/sfr_100'.ljust(128)
        glob_attrs.varnames.T.view('S128')[3] = b'star/metallicity'.ljust(128)
        # TODO: Add angular momentum computation
        # glob_attrs.varnames.T.view('S128')[3] = b'star/ang_momentum_x'.ljust(128)
        # glob_attrs.varnames.T.view('S128')[3] = b'star/ang_momentum_y'.ljust(128)
        # glob_attrs.varnames.T.view('S128')[3] = b'star/ang_momentum_z'.ljust(128)
        # glob_attrs.varnames.T.view('S128')[3] = b'dm/ang_momentum_x'.ljust(128)
        # glob_attrs.varnames.T.view('S128')[3] = b'dm/ang_momentum_y'.ljust(128)
        # glob_attrs.varnames.T.view('S128')[3] = b'dm/ang_momentum_z'.ljust(128)

        glob_attrs.varnames.T.view('S128')[4] = b'dm/mass'.ljust(128)
        glob_attrs.wvarnames.T.view('S128')[0] = b'cumulative'.ljust(128)
        glob_attrs.wvarnames.T.view('S128')[1] = b'mass'.ljust(128)

        # Begin integration
        part_integrator.integrate_region(output_path,selected_reg,filt,glob_attrs)

        print(glob_attrs.data)

        # DM details
        self.mass['dm'] = self.obj.yt_dataset.quan(glob_attrs.data[-1,0,0], 'code_mass')
        self.ndm = glob_attrs.ndm

        # Stellar details
        self.mass['stellar'] = self.obj.yt_dataset.quan(glob_attrs.data[0,0,0], 'code_mass')
        self.nstar = glob_attrs.nstar
        self.sfr['10Myr'] = self.obj.yt_dataset.quan(glob_attrs.data[1,0,0],'Msun/yr')
        self.sfr['100Myr'] = self.obj.yt_dataset.quan(glob_attrs.data[2,0,0],'Msun/yr')
        self.metallicity['stellar'] = glob_attrs.data[3,1,0] # Mass-weighted average!

    def _calculate_gas_quantities(self):
        """Compute gas quantities: Metallicity, Temperature..."""
        output_path = self.obj.simulation.fullpath

        # Initialise region
        selected_reg = init_region(self,'sphere')

        # We do not want any particular filter, just simple integration will do
        filt = filtering.filter()

        # Since the number of global quanties that can be computed
        # from the gas data depends on the specific configuration
        # of a simulation, this needs to be determined at the start
        # See: ozy/sim_attributes.py/assign_attributes

        nvar = 0
        quantity_names = []
        weight_names = ['cumulative','mass','volume']
        if self.obj.simulation.physics['hydro']:
            quantity_names += ['mass','density','temperature',
                                'ang_momentum_x','ang_momentum_y',
                                'ang_momentum_z','thermal_energy','thermal_energy_specific']
        if self.obj.simulation.physics['metals']:
            quantity_names += ['metallicity']
        if self.obj.simulation.physics['magnetic']:
            quantity_names += ['magnetic_energy','magnetic_energy_specific']
        if self.obj.simulation.physics['cr']:
            quantity_names += ['cr_energy','cr_energy_specific']
            
        nvar = len(quantity_names)
        if self.obj.simulation.physics['rt']:
            quantity_names += ['xHII','xHeII','xHeII']
                
        # Initialise Fortran derived type with attributes
        # This object hold the following attributes:
        # - nvars: number of variables
        # - nwvars:  number of variables for weighting
        # - varnames: names of variables
        # - wvarnames:  names of variables for weighting
        # - data: organised in numpy array of shape (nvars,nwvars,4)
        #           each of those 4 values are (final, min, max, sum of weights)
        glob_attrs = amr_integrator.amr_region_attrs()
        glob_attrs.nvars = len(quantity_names)
        glob_attrs.nwvars = len(weight_names)
        amr_integrator.allocate_amr_regions_attrs(glob_attrs)
        for i in range(0, len(quantity_names)):
            glob_attrs.varnames.T.view('S128')[i] = quantity_names[i].ljust(128)
        for i in range(0, len(weight_names)):
            glob_attrs.wvarnames.T.view('S128')[i] = weight_names[i].ljust(128)
        
        # Begin integration
        amr_integrator.integrate_region(output_path,selected_reg,filt,glob_attrs)

        # Assign results to galaxy object
        print(self.obj.simulation.physics)
        if self.obj.simulation.physics['hydro']:
            self.mass['gas'] = self.obj.yt_dataset.quan(glob_attrs.data[0,0,0], 'code_mass')
            self.gas_density['mass_weighted'] = self.obj.yt_dataset.quan(glob_attrs.data[1,1,0], 'code_density')
            self.gas_density['volume_weighted'] = self.obj.yt_dataset.quan(glob_attrs.data[1,2,0], 'code_density')
            self.temperature['mass_weighted'] = self.obj.yt_dataset.quan(glob_attrs.data[2,1,0], 'code_temperature')
            self.temperature['volume_weighted'] = self.obj.yt_dataset.quan(glob_attrs.data[2,2,0], 'code_temperature')
            self.angular_mom['gas'] = self.obj.yt_dataset.arr(np.array([glob_attrs.data[3,0,0],glob_attrs.data[4,0,0],glob_attrs.data[5,0,0]]),
                                                             'code_mass*code_length*code_velocity')
            self.energies['thermal_energy'] = self.obj.yt_dataset.quan(glob_attrs.data[6,0,0], 'code_mass * code_velocity**2')
            self.energies['thermal_energy_specific'] = self.obj.yt_dataset.quan(glob_attrs.data[7,2,0], 'code_specific_energy')

            if self.obj.simulation.physics['metals']:
                self.metallicity['gas'] = glob_attrs.data[8,1,0] # Mass-weighted average!
        else:
            self.mass['gas'] = self.obj.yt_dataset.quan(0.0, 'code_mass')
        
        if self.obj.simulation.physics['magnetic']:
            print('Computing magnetic energies')
            if self.obj.simulation.physics['metals']:
                self.energies['magnetic_energy'] = self.obj.yt_dataset.quan(glob_attrs.data[9,0,0], 'code_mass * code_velocity**2')
                self.energies['magnetic_energy_specific'] = self.obj.yt_dataset.quan(glob_attrs.data[10,1,0], 'code_specific_energy')
            else:
                self.energies['magnetic_energy'] = self.obj.yt_dataset.quan(glob_attrs.data[8,0,0], 'code_mass * code_velocity**2')
                self.energies['magnetic_energy_specific'] = self.obj.yt_dataset.quan(glob_attrs.data[9,1,0], 'code_specific_energy')

        if self.obj.simulation.physics['cr']:
            print('Computing CR energies')
            if self.obj.simulation.physics['metals']:
                self.energies['cr_energy'] = self.obj.yt_dataset.quan(glob_attrs.data[11,0,0], 'code_mass * code_velocity**2')
                self.energies['cr_energy_specific'] = self.obj.yt_dataset.quan(glob_attrs.data[12,1,0], 'code_specific_energy')
            else:
                self.energies['cr_energy'] = self.obj.yt_dataset.quan(glob_attrs.data[10,0,0], 'code_mass * code_velocity**2')
                self.energies['cr_energy_specific'] = self.obj.yt_dataset.quan(glob_attrs.data[11,1,0], 'code_specific_energy')

        if self.obj.simulation.physics['rt']:
            print('Computing ionisation fractions')
            self.radiation['xHII'] = glob_attrs.data[nvar,1,0]
            self.radiation['xHeII'] = glob_attrs.data[nvar+1,1,0]
            self.radiation['xHeIII'] = glob_attrs.data[nvar+2,1,0]
            
    def _calculate_star_quantities(self):
        """Calculate star quantities..."""
        # TODO
        indicators = np.array([0.01, 0.1]) # in Gyr. These should be taken as arguments in the future
        output_path = self.obj.simulation.fullpath

        # Region for selection -- Using: 0.2 radius of host DM halo (following Martin-Alvarez et al. 2018)
        r_region = 0.2*self.obj.halos[self.parent_halo_index].virial_quantities['radius'].in_units('code_length').d
        x_center = self.position.in_units('code_length')[0].d
        y_center = self.position.in_units('code_length')[1].d
        z_center = self.position.in_units('code_length')[2].d
        
        sfr = part2sfr.sphere(output_path, x_center, y_center,
                                z_center, r_region, indicators, 
                                len(indicators))

        for i in range(0, len(indicators)):
            self.sfr[str(int(indicators[i]*1e+3))+'Myr'] = self.obj.yt_dataset.quan(sfr[i],'Msun/yr')

    def _calculate_velocity_dispersions(self):
        """Calculate velocity dispersions for the various components."""
        # TODO
        return
class Halo(Group):
    """Halo class which has different levels of the halo hierarchy."""
    obj_type = 'halo'
    def __init__(self, obj):
        super(Halo, self).__init__(obj)
        self.spin = 0
        self.level              = -1
        self.host               = -1
        self.hostsub            = -1
        self.nsub               = -1
        self.nextsub            = -1
        self.galaxies           = []
        self.central_galaxy     = None
        self.satellite_galaxies = []
        self.galaxy_index_list = []

def create_new_group(obj, grouptype):
    """Simple function to create a new instance of a specified :class:`group.Group`.
    """
    if grouptype == 'halo':
        return Halo(obj)
    elif grouptype == 'galaxy':
        return Galaxy(obj)