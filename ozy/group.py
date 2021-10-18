import numpy as np
from pprint import pprint
import sys
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/amr')
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/part')
from amr2 import vectors
from amr2 import geometrical_regions as geo
from amr2 import filtering
from amr2 import amr_integrator
from part2 import part_integrator

from ozy.profiles import init_region,init_filter

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
        self.outflows = {}
        self.inflows = {}
    def _process_galaxy(self):
        """Process each galaxy after creation. This means
        calculating the total mass, and then calculate the rest of masses,
        radial quantities, velocity dispersions, angular momentum...
        """
        self._calculate_stardm_quantities()
        self._calculate_velocity_dispersions()
        self._calculate_gas_quantities()
        # self._calculate_outflow_inflow()

        self.mass['baryon'] = self.mass['stellar'] + self.mass['gas']

        self.gas_fraction = 0.0
        if self.mass['baryon'] > 0:
            self.gas_fraction = self.mass['gas'].d / self.mass['baryon'].d
    def _empty_galaxy(self):
        """Add atributes to galaxy to zero when no interest is in it."""

        self.mass['dm'] = self.obj.yt_dataset.quan(0.0, 'code_mass')
        self.mass['stellar'] = self.obj.yt_dataset.quan(0.0, 'code_mass')
        self.mass['gas'] = self.obj.yt_dataset.quan(0.0, 'code_mass')
        self.mass['baryon'] = self.mass['stellar'] + self.mass['gas']
        self.gas_fraction = 0.0

        self.ndm = 0

        # Stellar details
        self.mass['stellar'] = self.obj.yt_dataset.quan(0.0, 'code_mass')
        self.nstar = 0
        self.sfr['10Myr'] = self.obj.yt_dataset.quan(0.0,'Msun/yr')
        self.sfr['100Myr'] = self.obj.yt_dataset.quan(0.0,'Msun/yr')
        self.metallicity['stellar'] = 0 # Mass-weighted average!

        if self.obj.simulation.physics['hydro']:
            self.mass['gas'] = self.obj.yt_dataset.quan(0.0, 'code_mass')
            self.gas_density['mass_weighted'] = self.obj.yt_dataset.quan(0.0, 'code_density')
            self.gas_density['volume_weighted'] = self.obj.yt_dataset.quan(0.0, 'code_density')
            self.temperature['mass_weighted'] = self.obj.yt_dataset.quan(0.0, 'code_temperature')
            self.temperature['volume_weighted'] = self.obj.yt_dataset.quan(0.0, 'code_temperature')
            self.angular_mom['gas'] = self.obj.yt_dataset.arr(np.array([0,0,0]),
                                                             'code_mass*code_length*code_velocity')
            self.energies['thermal_energy'] = self.obj.yt_dataset.quan(0, 'code_mass * code_velocity**2')
            self.energies['thermal_energy_specific'] = self.obj.yt_dataset.quan(0, 'code_specific_energy')

            if self.obj.simulation.physics['metals']:
                self.metallicity['gas'] = 0 # Mass-weighted average!
        else:
            self.mass['gas'] = self.obj.yt_dataset.quan(0.0, 'code_mass')
        
        if self.obj.simulation.physics['magnetic']:
            if self.obj.simulation.physics['metals']:
                self.energies['magnetic_energy'] = self.obj.yt_dataset.quan(0, 'code_mass * code_velocity**2')
                self.energies['magnetic_energy_specific'] = self.obj.yt_dataset.quan(0, 'code_specific_energy')
            else:
                self.energies['magnetic_energy'] = self.obj.yt_dataset.quan(0, 'code_mass * code_velocity**2')
                self.energies['magnetic_energy_specific'] = self.obj.yt_dataset.quan(0, 'code_specific_energy')

        if self.obj.simulation.physics['cr']:
            if self.obj.simulation.physics['metals']:
                self.energies['cr_energy'] = self.obj.yt_dataset.quan(0, 'code_mass * code_velocity**2')
                self.energies['cr_energy_specific'] = self.obj.yt_dataset.quan(0, 'code_specific_energy')
            else:
                self.energies['cr_energy'] = self.obj.yt_dataset.quan(0, 'code_mass * code_velocity**2')
                self.energies['cr_energy_specific'] = self.obj.yt_dataset.quan(0, 'code_specific_energy')

        if self.obj.simulation.physics['rt']:
            self.radiation['xHII'] = 0
            self.radiation['xHeII'] = 0
            self.radiation['xHeIII'] = 0


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
        glob_attrs.nvars = 11
        glob_attrs.nwvars = 2
        part_integrator.allocate_part_regions_attrs(glob_attrs)
        glob_attrs.varnames.T.view('S128')[0] = b'star/mass'.ljust(128)
        # TODO: SFR indicators should be taken as arguments
        glob_attrs.varnames.T.view('S128')[1] = b'star/sfr_10'.ljust(128)
        glob_attrs.varnames.T.view('S128')[2] = b'star/sfr_100'.ljust(128)
        glob_attrs.varnames.T.view('S128')[3] = b'star/metallicity'.ljust(128)

        glob_attrs.varnames.T.view('S128')[4] = b'star/ang_momentum_x'.ljust(128)
        glob_attrs.varnames.T.view('S128')[5] = b'star/ang_momentum_y'.ljust(128)
        glob_attrs.varnames.T.view('S128')[6] = b'star/ang_momentum_z'.ljust(128)

        glob_attrs.varnames.T.view('S128')[7] = b'dm/mass'.ljust(128)
        glob_attrs.varnames.T.view('S128')[8] = b'dm/ang_momentum_x'.ljust(128)
        glob_attrs.varnames.T.view('S128')[9] = b'dm/ang_momentum_y'.ljust(128)
        glob_attrs.varnames.T.view('S128')[10] = b'dm/ang_momentum_z'.ljust(128)

        
        glob_attrs.wvarnames.T.view('S128')[0] = b'cumulative'.ljust(128)
        glob_attrs.wvarnames.T.view('S128')[1] = b'mass'.ljust(128)

        # Begin integration
        part_integrator.integrate_region(output_path,selected_reg,filt,glob_attrs)

        # Stellar details
        self.mass['stellar'] = self.obj.yt_dataset.quan(glob_attrs.data[0,0,0], 'code_mass')
        self.nstar = glob_attrs.nstar
        self.sfr['10Myr'] = self.obj.yt_dataset.quan(glob_attrs.data[1,0,0],'Msun/yr')
        self.sfr['100Myr'] = self.obj.yt_dataset.quan(glob_attrs.data[2,0,0],'Msun/yr')
        self.metallicity['stellar'] = glob_attrs.data[3,1,0] # Mass-weighted average!
        self.angular_mom['stellar'] = self.obj.yt_dataset.arr(np.array([glob_attrs.data[4,0,0],glob_attrs.data[5,0,0],glob_attrs.data[6,0,0]]),
                                                             'code_mass*code_length*code_velocity')

        # DM details
        self.mass['dm'] = self.obj.yt_dataset.quan(glob_attrs.data[7,0,0], 'code_mass')
        self.ndm = glob_attrs.ndm
        self.angular_mom['dm'] = self.obj.yt_dataset.arr(np.array([glob_attrs.data[8,0,0],glob_attrs.data[9,0,0],glob_attrs.data[10,0,0]]),
                                                             'code_mass*code_length*code_velocity')


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
            quantity_names += ['xHII','xHeII','xHeIII']
                
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
            self.radiation['xHII'] = glob_attrs.data[nvar,2,0]
            self.radiation['xHeII'] = glob_attrs.data[nvar+1,2,0]
            self.radiation['xHeIII'] = glob_attrs.data[nvar+2,2,0]

    def _calculate_outflow_inflow(self):
        """Compute details of outflows and inflows, by measuring quantities on a thin shell."""
        output_path = self.obj.simulation.fullpath

        # Define shell regions at different radii with width of 0.02*r_vir
        shell_width = 0.02*self.obj.halos[self.parent_halo_index].virial_quantities['radius'].d
        shell_02 = init_region(self, 'sphere',rmin=0.19,rmax=0.21)
        shell_05 = init_region(self, 'sphere',rmin=0.49,rmax=0.51)
        shell_10 = init_region(self, 'sphere',rmin=0.99,rmax=1.01)

        # Initialise filter for inflow and outflow gas
        out_filter = init_filter(cond_strs='v_sphere_r/>/0/km*s**-1',name='outflows',group=self)
        in_filter = init_filter(cond_strs='v_sphere_r/<=/0/km*s**-1',name='inflows',group=self)

        # Since the number of global quanties that can be computed
        # from the gas data depends on the specific configuration
        # of a simulation, this needs to be determined at the start
        # See: ozy/sim_attributes.py/assign_attributes

        nvar = 0
        quantity_names = []
        weight_names = ['cumulative','massflux_rate_sphere_r']

        if self.obj.simulation.physics['hydro']:
            quantity_names += ['density','temperature',
                                'momentum_sphere_r','v_sphere_r',
                                'thermal_energy','thermal_energy_specific']
        if self.obj.simulation.physics['metals']:
            quantity_names += ['metallicity']
        if self.obj.simulation.physics['magnetic']:
            quantity_names += ['magnetic_energy','magnetic_energy_specific']
        if self.obj.simulation.physics['cr']:
            quantity_names += ['cr_energy','cr_energy_specific']

        nvar = len(quantity_names)
        if self.obj.simulation.physics['rt']:
            quantity_names += ['xHII','xHeII','xHeIII']

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

        # First shell
        print('Performing first shell integration')
        # a) Outflows
        amr_integrator.integrate_region(output_path,shell_02,out_filter,glob_attrs)

        # Assign results to galaxy object
        if self.obj.simulation.physics['hydro']:
            print('Computing gas flow quantities')
            self.outflows['density_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[0,1,0], 'code_density')
            self.outflows['temperature_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[1,1,0], 'code_temperature')
            self.outflows['massflow_rate_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[2,0,0]/shell_width, 'code_mass*code_velocity/code_length')
            self.outflows['v_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[3,1,0], 'code_velocity')
            self.outflows['thermal_energy_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[4,0,0], 'code_mass * code_velocity**2')
            self.outflows['thermal_energy_specific_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[5,1,0], 'code_specific_energy')
            print(self.outflows)
            if self.obj.simulation.physics['metals']:
                self.outflows['metallicity_20rvir'] = glob_attrs.data[6,1,0]
        else:
            self.outflows['massflow_rate_20rvir'] = self.obj.yt_dataset.quan(0.0, 'code_mass*code_velocity/code_length')
        
        if self.obj.simulation.physics['magnetic']:
            print('Computing magnetic energies')
            if self.obj.simulation.physics['metals']:
                self.outflows['magnetic_energy_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[7,0,0], 'code_mass * code_velocity**2')
                self.outflows['magnetic_energy_specific_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[8,1,0], 'code_specific_energy')
            else:
                self.outflows['magnetic_energy_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[6,0,0], 'code_mass * code_velocity**2')
                self.outflows['magnetic_energy_specific_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[7,1,0], 'code_specific_energy')

        if self.obj.simulation.physics['cr']:
            print('Computing CR energies')
            if self.obj.simulation.physics['metals']:
                self.outflows['cr_energy_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[9,0,0], 'code_mass * code_velocity**2')
                self.outflows['cr_energy_specific_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[10,1,0], 'code_specific_energy')
            else:
                self.outflows['cr_energy_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[8,0,0], 'code_mass * code_velocity**2')
                self.outflows['cr_energy_specific_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[9,1,0], 'code_specific_energy')

        if self.obj.simulation.physics['rt']:
            print('Computing ionisation fractions')
            self.outflows['xHII_20rvir'] = glob_attrs.data[nvar,2,0]
            self.outflows['xHeII_20rvir'] = glob_attrs.data[nvar+1,2,0]
            self.outflows['xHeIII_20rvir'] = glob_attrs.data[nvar+2,2,0]

        # b) Inflows
        amr_integrator.integrate_region(output_path,shell_02,in_filter,glob_attrs)

        # Assign results to galaxy object
        if self.obj.simulation.physics['hydro']:
            print('Computing gas flow quantities')
            self.inflows['density_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[0,1,0], 'code_density')
            self.inflows['temperature_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[1,1,0], 'code_temperature')
            self.inflows['massflow_rate_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[2,0,0]/shell_width, 'code_mass*code_velocity/code_length')
            self.inflows['v_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[3,1,0], 'code_velocity')
            self.inflows['thermal_energy_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[4,0,0], 'code_mass * code_velocity**2')
            self.inflows['thermal_energy_specific_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[5,1,0], 'code_specific_energy')

            if self.obj.simulation.physics['metals']:
                self.outflows['metallicity_20rvir'] = glob_attrs.data[6,1,0]
        else:
            self.inflows['massflow_rate_20rvir'] = self.obj.yt_dataset.quan(0.0, 'code_mass*code_velocity/code_length')
        
        if self.obj.simulation.physics['magnetic']:
            print('Computing magnetic energies')
            if self.obj.simulation.physics['metals']:
                self.inflows['magnetic_energy_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[7,0,0], 'code_mass * code_velocity**2')
                self.inflows['magnetic_energy_specific_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[8,1,0], 'code_specific_energy')
            else:
                self.inflows['magnetic_energy_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[6,0,0], 'code_mass * code_velocity**2')
                self.inflows['magnetic_energy_specific_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[7,1,0], 'code_specific_energy')

        if self.obj.simulation.physics['cr']:
            print('Computing CR energies')
            if self.obj.simulation.physics['metals']:
                self.inflows['cr_energy_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[9,0,0], 'code_mass * code_velocity**2')
                self.inflows['cr_energy_specific_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[10,1,0], 'code_specific_energy')
            else:
                self.inflows['cr_energy_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[8,0,0], 'code_mass * code_velocity**2')
                self.inflows['cr_energy_specific_20rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[9,1,0], 'code_specific_energy')

        if self.obj.simulation.physics['rt']:
            print('Computing ionisation fractions')
            self.inflows['xHII_20rvir'] = glob_attrs.data[nvar,2,0]
            self.inflows['xHeII_20rvir'] = glob_attrs.data[nvar+1,2,0]
            self.inflows['xHeIII_20rvir'] = glob_attrs.data[nvar+2,2,0]

        # Second shell
        print('Performing second shell integration')
        # a) Outflows
        amr_integrator.integrate_region(output_path,shell_05,out_filter,glob_attrs)

        # Assign results to galaxy object
        if self.obj.simulation.physics['hydro']:
            print('Computing gas flow quantities')
            self.outflows['density_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[0,1,0], 'code_density')
            self.outflows['temperature_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[1,1,0], 'code_temperature')
            self.outflows['massflow_rate_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[2,0,0]/shell_width, 'code_mass*code_velocity/code_length')
            self.outflows['v_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[3,1,0], 'code_velocity')
            self.outflows['thermal_energy_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[4,0,0], 'code_mass * code_velocity**2')
            self.outflows['thermal_energy_specific_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[5,1,0], 'code_specific_energy')

            if self.obj.simulation.physics['metals']:
                self.outflows['metallicity_50rvir'] = glob_attrs.data[6,1,0]
        else:
            self.outflows['massflow_rate_50rvir'] = self.obj.yt_dataset.quan(0.0, 'code_mass*code_velocity/code_length')
        
        if self.obj.simulation.physics['magnetic']:
            print('Computing magnetic energies')
            if self.obj.simulation.physics['metals']:
                self.outflows['magnetic_energy_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[7,0,0], 'code_mass * code_velocity**2')
                self.outflows['magnetic_energy_specific_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[8,1,0], 'code_specific_energy')
            else:
                self.outflows['magnetic_energy_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[6,0,0], 'code_mass * code_velocity**2')
                self.outflows['magnetic_energy_specific_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[7,1,0], 'code_specific_energy')

        if self.obj.simulation.physics['cr']:
            print('Computing CR energies')
            if self.obj.simulation.physics['metals']:
                self.outflows['cr_energy_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[9,0,0], 'code_mass * code_velocity**2')
                self.outflows['cr_energy_specific_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[10,1,0], 'code_specific_energy')
            else:
                self.outflows['cr_energy_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[8,0,0], 'code_mass * code_velocity**2')
                self.outflows['cr_energy_specific_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[9,1,0], 'code_specific_energy')

        if self.obj.simulation.physics['rt']:
            print('Computing ionisation fractions')
            self.outflows['xHII_50rvir'] = glob_attrs.data[nvar,2,0]
            self.outflows['xHeII_50rvir'] = glob_attrs.data[nvar+1,2,0]
            self.outflows['xHeIII_50rvir'] = glob_attrs.data[nvar+2,2,0]

        # b) Inflows
        amr_integrator.integrate_region(output_path,shell_05,in_filter,glob_attrs)

        # Assign results to galaxy object
        if self.obj.simulation.physics['hydro']:
            print('Computing gas flow quantities')
            self.inflows['density_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[0,1,0], 'code_density')
            self.inflows['temperature_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[1,1,0], 'code_temperature')
            self.inflows['massflow_rate_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[2,0,0]/shell_width, 'code_mass*code_velocity/code_length')
            self.inflows['v_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[3,1,0], 'code_velocity')
            self.inflows['thermal_energy_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[4,0,0], 'code_mass * code_velocity**2')
            self.inflows['thermal_energy_specific_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[5,1,0], 'code_specific_energy')

            if self.obj.simulation.physics['metals']:
                self.outflows['metallicity_50rvir'] = glob_attrs.data[6,1,0]
        else:
            self.inflows['massflow_rate_50rvir'] = self.obj.yt_dataset.quan(0.0, 'code_mass*code_velocity/code_length')
        
        if self.obj.simulation.physics['magnetic']:
            print('Computing magnetic energies')
            if self.obj.simulation.physics['metals']:
                self.inflows['magnetic_energy_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[7,0,0], 'code_mass * code_velocity**2')
                self.inflows['magnetic_energy_specific_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[8,1,0], 'code_specific_energy')
            else:
                self.inflows['magnetic_energy_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[6,0,0], 'code_mass * code_velocity**2')
                self.inflows['magnetic_energy_specific_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[7,1,0], 'code_specific_energy')

        if self.obj.simulation.physics['cr']:
            print('Computing CR energies')
            if self.obj.simulation.physics['metals']:
                self.inflows['cr_energy_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[9,0,0], 'code_mass * code_velocity**2')
                self.inflows['cr_energy_specific_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[10,1,0], 'code_specific_energy')
            else:
                self.inflows['cr_energy_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[8,0,0], 'code_mass * code_velocity**2')
                self.inflows['cr_energy_specific_50rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[9,1,0], 'code_specific_energy')

        if self.obj.simulation.physics['rt']:
            print('Computing ionisation fractions')
            self.inflows['xHII_50rvir'] = glob_attrs.data[nvar,2,0]
            self.inflows['xHeII_50rvir'] = glob_attrs.data[nvar+1,2,0]
            self.inflows['xHeIII_50rvir'] = glob_attrs.data[nvar+2,2,0]

        # Third shell
        print('Performing third shell integration')
        # a) Outflows
        amr_integrator.integrate_region(output_path,shell_10,out_filter,glob_attrs)

        # Assign results to galaxy object
        if self.obj.simulation.physics['hydro']:
            print('Computing gas flow quantities')
            self.outflows['density_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[0,1,0], 'code_density')
            self.outflows['temperature_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[1,1,0], 'code_temperature')
            self.outflows['massflow_rate_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[2,0,0]/shell_width, 'code_mass*code_velocity/code_length')
            self.outflows['v_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[3,1,0], 'code_velocity')
            self.outflows['thermal_energy_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[4,0,0], 'code_mass * code_velocity**2')
            self.outflows['thermal_energy_specific_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[5,1,0], 'code_specific_energy')

            if self.obj.simulation.physics['metals']:
                self.outflows['metallicity_100rvir'] = glob_attrs.data[6,1,0]
        else:
            self.outflows['massflow_rate_100rvir'] = self.obj.yt_dataset.quan(0.0, 'code_mass*code_velocity/code_length')
        
        if self.obj.simulation.physics['magnetic']:
            print('Computing magnetic energies')
            if self.obj.simulation.physics['metals']:
                self.outflows['magnetic_energy_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[7,0,0], 'code_mass * code_velocity**2')
                self.outflows['magnetic_energy_specific_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[8,1,0], 'code_specific_energy')
            else:
                self.outflows['magnetic_energy_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[6,0,0], 'code_mass * code_velocity**2')
                self.outflows['magnetic_energy_specific_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[7,1,0], 'code_specific_energy')

        if self.obj.simulation.physics['cr']:
            print('Computing CR energies')
            if self.obj.simulation.physics['metals']:
                self.outflows['cr_energy_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[9,0,0], 'code_mass * code_velocity**2')
                self.outflows['cr_energy_specific_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[10,1,0], 'code_specific_energy')
            else:
                self.outflows['cr_energy_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[8,0,0], 'code_mass * code_velocity**2')
                self.outflows['cr_energy_specific_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[9,1,0], 'code_specific_energy')

        if self.obj.simulation.physics['rt']:
            print('Computing ionisation fractions')
            self.outflows['xHII_100rvir'] = glob_attrs.data[nvar,2,0]
            self.outflows['xHeII_100rvir'] = glob_attrs.data[nvar+1,2,0]
            self.outflows['xHeIII_100rvir'] = glob_attrs.data[nvar+2,2,0]

        # b) Inflows
        amr_integrator.integrate_region(output_path,shell_10,in_filter,glob_attrs)

        # Assign results to galaxy object
        if self.obj.simulation.physics['hydro']:
            print('Computing gas flow quantities')
            self.inflows['density_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[0,1,0], 'code_density')
            self.inflows['temperature_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[1,1,0], 'code_temperature')
            self.inflows['massflow_rate_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[2,0,0]/shell_width, 'code_mass*code_velocity/code_length')
            self.inflows['v_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[3,1,0], 'code_velocity')
            self.inflows['thermal_energy_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[4,0,0], 'code_mass * code_velocity**2')
            self.inflows['thermal_energy_specific_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[5,1,0], 'code_specific_energy')

            if self.obj.simulation.physics['metals']:
                self.outflows['metallicity_100rvir'] = glob_attrs.data[6,1,0]
        else:
            self.inflows['massflow_rate_100rvir'] = self.obj.yt_dataset.quan(0.0, 'code_mass*code_velocity/code_length')
        
        if self.obj.simulation.physics['magnetic']:
            print('Computing magnetic energies')
            if self.obj.simulation.physics['metals']:
                self.inflows['magnetic_energy_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[7,0,0], 'code_mass * code_velocity**2')
                self.inflows['magnetic_energy_specific_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[8,1,0], 'code_specific_energy')
            else:
                self.inflows['magnetic_energy_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[6,0,0], 'code_mass * code_velocity**2')
                self.inflows['magnetic_energy_specific_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[7,1,0], 'code_specific_energy')

        if self.obj.simulation.physics['cr']:
            print('Computing CR energies')
            if self.obj.simulation.physics['metals']:
                self.inflows['cr_energy_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[9,0,0], 'code_mass * code_velocity**2')
                self.inflows['cr_energy_specific_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[10,1,0], 'code_specific_energy')
            else:
                self.inflows['cr_energy_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[8,0,0], 'code_mass * code_velocity**2')
                self.inflows['cr_energy_specific_100rvir'] = self.obj.yt_dataset.quan(glob_attrs.data[9,1,0], 'code_specific_energy')

        if self.obj.simulation.physics['rt']:
            print('Computing ionisation fractions')
            self.inflows['xHII_100rvir'] = glob_attrs.data[nvar,2,0]
            self.inflows['xHeII_100rvir'] = glob_attrs.data[nvar+1,2,0]
            self.inflows['xHeIII_100rvir'] = glob_attrs.data[nvar+2,2,0]
        
        print('Outflow/inflow analysis finished!')
        return

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
