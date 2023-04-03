import numpy as np
from pprint import pprint
from amr2 import vectors
from amr2 import geometrical_regions as geo
from amr2 import filtering
from amr2 import amr_integrator,stats_utils
from part2 import part_integrator

from ozy.utils import init_region,init_filter

MINIMUM_STARS_PER_GALAXY = 100
MINIMUM_DM_PER_HALO      = 0

grouptypes = dict(
    halo='halos',
    galaxy='galaxies',
    region='regions',
)

info_blacklist = [
    'obj', 'halo', 'galaxies', 'satellites',
    'galaxy_index_list_end', 'galaxy_index_list_start']

class Group(object):
    """This is the parent class from which the rest of groups are derived."""
    
    def __init__(self, obj):
        """Initialise basic properties of group object."""
        self.obj      = obj
        self.ID       = -1
        self.level    = -1
        self.host     = -1
        self.hostsub  = -1
        self.nsub     = -1
        self.nextsub  = -1
        self.npart    = 0
        self.units    = 'code'
        self.position = np.array([-1.0, -1.0, -1.0])
        self.radius   = {}
        self.shape = {}
        self.mass     = {}
        self.velocity = np.array([0.0, 0.0, 0.0])
        self.angular_mom = {}
        self.virial_quantities = {}
        self.energies = {}

        if 'tidal_method' in self.obj._kwargs:
            tidal_method  = self.obj._kwargs['tidal_method']
        else:
            tidal_method = 'BT87_simple'
        self.radius[tidal_method] = self.obj.quantity(0.0, 'code_length')
    
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
        self.magnetism = {}
        self.sf_efficiency = {}
        self.pressure_support = {}
        self.velocity_dispersion = {}
    def _process_galaxy(self):
        """Process each galaxy after creation. This means
        calculating the total mass, and then calculate the rest of masses,
        radial quantities, velocity dispersions, angular momentum...
        """
        self._calculate_stardm_quantities()
        self._calculate_velocity_dispersions()
        self._calculate_galaxy_gas_quantities()
        self._calculate_halo_gas_quantities()
        if 'tidal_method' in self.obj._kwargs:
            tidal_method  = self.obj._kwargs['tidal_method']
        else:
            tidal_method = 'BT87_simple'
        self._calculate_substructure_tidal(tidal_method)

        self.mass['baryon'] = self.mass['stellar'] + self.mass['gas']

        self.gas_fraction = 0.0
        if self.mass['baryon'] > 0:
            self.gas_fraction = self.mass['gas'].d / self.mass['baryon'].d
    def _empty_galaxy(self):
        """Add atributes to galaxy to zero when no interest is in it."""

        self.mass['dm'] = self.obj.quantity(0.0, 'code_mass')
        self.mass['stellar'] = self.obj.quantity(0.0, 'code_mass')
        self.mass['gas'] = self.obj.quantity(0.0, 'code_mass')
        self.mass['baryon'] = self.mass['stellar'] + self.mass['gas']
        self.gas_fraction = 0.0

        self.ndm = 0

        # Stellar details
        self.mass['stellar'] = self.obj.quantity(0.0, 'code_mass')
        self.nstar = 0
        self.sfr['10Myr'] = self.obj.quantity(0.0,'Msun/yr')
        self.sfr['100Myr'] = self.obj.quantity(0.0,'Msun/yr')
        self.metallicity['stellar'] = 0 # Mass-weighted average!

        phase_names = ['cold','warm','hot']

        empty_array = np.zeros((4,7)).astype(np.float64)

        if self.obj.simulation.physics['hydro']:
            self.mass['gas'] = self.obj.quantity(0.0, 'code_mass')
            self.gas_density['galaxy_gas'] = self.obj.array(empty_array, 'code_density')
            self.temperature['galaxy_gas'] = self.obj.array(empty_array, 'code_temperature')
            self.angular_mom['gas'] = self.obj.array(np.array([0.0,0.0,0.0]).astype(np.float64),
                                                             'code_mass*code_length*code_velocity')
            self.energies['thermal_energy'] = self.obj.quantity(0.0, 'code_mass * code_velocity**2')
            self.energies['thermal_energy_specific'] = self.obj.array(empty_array, 'code_specific_energy')
            self.pressure_support['grav_therpfrsphere'] = self.obj.array(empty_array, 'dimensionless')
            self.velocity_dispersion['gas_turbulent'] = self.obj.array(empty_array, 'dimensionless')
            if not self.obj.simulation.physics['magnetic'] and not self.obj.simulation.physics['cr']:
                self.sf_efficiency['eff_FK2'] = self.obj.array(empty_array, 'dimensionless')
            if self.obj.simulation.physics['metals']:
                self.metallicity['galaxy_gas'] = self.obj.array(empty_array, 'dimensionless')
            for i in range(0, len(phase_names)):
                self.mass['gas_'+phase_names[i]] = self.obj.quantity(0.0, 'code_mass')
                self.gas_density[phase_names[i]] = self.obj.array(empty_array, 'code_density')
                self.temperature[phase_names[i]] = self.obj.array(empty_array, 'code_temperature')
                self.angular_mom['gas_'+phase_names[i]] = self.obj.array(np.array([0.0,0.0,0.0]).astype(np.float64),
                                                                        'code_mass*code_length*code_velocity')
                self.energies['thermal_energy_'+phase_names[i]] = self.obj.quantity(0.0, 'code_mass * code_velocity**2')
                self.energies['thermal_energy_specific_'+phase_names[i]] = self.obj.array(empty_array, 'code_specific_energy')
                self.pressure_support['grav_therpfrsphere_'+phase_names[i]] = self.obj.array(empty_array, 'dimensionless')
                self.velocity_dispersion['gas_turbulent_'+phase_names[i]] = self.obj.array(empty_array, 'code_velocity')
                if self.obj.simulation.physics['metals']:
                    self.metallicity['gas_'+phase_names[i]] = self.obj.array(empty_array,'dimensionless')
        else:
            self.mass['gas'] = self.obj.quantity(0.0, 'code_mass')
        
        if self.obj.simulation.physics['magnetic']:
            self.energies['magnetic_energy'] = self.obj.quantity(0.0, 'code_mass * code_velocity**2')
            self.energies['magnetic_energy_specific'] = self.obj.array(empty_array, 'code_specific_energy')
            self.magnetism['magnetic_magnitude'] = self.obj.array(empty_array, 'code_magnetic')
            if not self.obj.simulation.physics['cr']:
                self.sf_efficiency['eff_FKmag'] = self.obj.array(empty_array, 'dimensionless')
            for i in range(0,len(phase_names)):
                self.energies['magnetic_energy_'+phase_names[i]] = self.obj.quantity(0.0, 'code_mass * code_velocity**2')
                self.energies['magnetic_energy_specific_'+phase_names[i]] = self.obj.array(empty_array, 'code_specific_energy')
                self.magnetism['magnetic_magnitude_'+phase_names[i]] = self.obj.array(empty_array, 'code_magnetic')

        if self.obj.simulation.physics['cr']:
            self.energies['cr_energy'] = self.obj.quantity(0.0, 'code_mass * code_velocity**2')
            self.energies['cr_energy_specific'] = self.obj.array(empty_array, 'code_specific_energy')
            self.pressure_support['grav_crpfrsphere'] = self.obj.array(empty_array, 'dimensionless')
            self.sf_efficiency['eff_FKmag'] = self.obj.array(empty_array, 'dimensionless')
            self.sf_efficiency['eff_FKmagnocr'] = self.obj.array(empty_array, 'dimensionless')
            for i in range(0, len(phase_names)):
                self.energies['cr_energy_'+phase_names[i]] = self.obj.quantity(0.0, 'code_mass * code_velocity**2')
                self.energies['cr_energy_specific_'+phase_names[i]] = self.obj.array(empty_array, 'code_specific_energy')
                self.pressure_support['grav_crpfrsphere_'+phase_names[i]] = self.obj.array(empty_array, 'dimensionless')
        if self.obj.simulation.physics['rt']:
            self.radiation['xHII'] = empty_array
            self.radiation['xHeII'] = empty_array
            self.radiation['xHeIII'] = empty_array
            for i in range(0, len(phase_names)):
                self.radiation['xHII_'+phase_names[i]] = empty_array
                self.radiation['xHeII_'+phase_names[i]] = empty_array
                self.radiation['xHeIII_'+phase_names[i]] = empty_array

        # Assign results to galaxy object
        if self.obj.simulation.physics['hydro']:
            self.mass['halo_gas'] = self.obj.quantity(0.0, 'code_mass')
            self.gas_density['halo_gas'] = self.obj.array(empty_array, 'code_density')
            self.temperature['halo_gas'] = self.obj.array(empty_array, 'code_temperature')
            self.angular_mom['halo_gas'] = self.obj.array(np.array([0.0,0.0,0.0]).astype(np.float64),
                                                             'code_mass*code_length*code_velocity')
            self.energies['halo_thermal_energy'] = self.obj.quantity(0.0, 'code_mass * code_velocity**2')
            self.energies['halo_thermal_energy_specific'] = self.obj.array(empty_array, 'code_specific_energy')
            self.pressure_support['halo_grav_therpfrsphere'] = self.obj.array(empty_array, 'dimensionless')
            self.velocity_dispersion['halo_gas_turbulent'] = self.obj.array(empty_array, 'code_velocity')
            if self.obj.simulation.physics['metals']:
                self.metallicity['halo_gas'] = self.obj.array(empty_array,'dimensionless')
        else:
            self.mass['halo_gas'] = self.obj.quantity(0.0, 'code_mass')
        
        if self.obj.simulation.physics['magnetic']:
            self.energies['halo_magnetic_energy'] = self.obj.quantity(0.0, 'code_mass * code_velocity**2')
            self.energies['halo_magnetic_energy_specific'] = self.obj.array(empty_array, 'code_specific_energy')
            self.magnetism['halo_magnetic_magnitude'] = self.obj.array(empty_array, 'code_magnetic')
        
        if self.obj.simulation.physics['cr']:
            self.energies['halo_cr_energy'] = self.obj.quantity(0.0, 'code_mass * code_velocity**2')
            self.energies['halo_cr_energy_specific'] = self.obj.array(empty_array, 'code_specific_energy')
            self.pressure_support['grav_crpfrsphere'] = self.obj.array(empty_array, 'dimensionless')

        if self.obj.simulation.physics['rt']:
            print('Computing ionisation fractions')
            self.radiation['halo_xHII'] = empty_array
            self.radiation['halo_xHeII'] = empty_array
            self.radiation['halo_xHeIII'] = empty_array


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
        self.mass['stellar'] = self.obj.quantity(glob_attrs.data[0,0,0], 'code_mass')
        self.nstar = glob_attrs.nstar
        self.sfr['10Myr'] = self.obj.quantity(glob_attrs.data[1,0,0],'Msun/yr')
        self.sfr['100Myr'] = self.obj.quantity(glob_attrs.data[2,0,0],'Msun/yr')
        self.metallicity['stellar'] = glob_attrs.data[3,1,0] # Mass-weighted average!
        self.angular_mom['stellar'] = self.obj.array(np.array([glob_attrs.data[4,0,0],glob_attrs.data[5,0,0],glob_attrs.data[6,0,0]]),
                                                             'code_mass*code_length*code_velocity')

        # DM details
        self.mass['dm'] = self.obj.quantity(glob_attrs.data[7,0,0], 'code_mass')
        self.ndm = glob_attrs.ndm
        self.angular_mom['dm'] = self.obj.array(np.array([glob_attrs.data[8,0,0],glob_attrs.data[9,0,0],glob_attrs.data[10,0,0]]),
                                                             'code_mass*code_length*code_velocity')


    def _calculate_galaxy_gas_quantities(self):
        from ozy.utils import get_code_bins, pdf_handler_to_stats
        """Compute gas quantities: Metallicity, Temperature..."""
        output_path = self.obj.simulation.fullpath

        # Initialise region
        selected_reg = init_region(self,'sphere',rmin=(0.0,'rvir'),rmax=(0.2,'rvir'))

        # We do not want any particular filter, just simple integration will do
        all_filt = filtering.filter()
        phase_names = ['cold','warm','hot']
        cold_filt = init_filter(cond_strs=['entropy_specific/</4.4e+8/erg*K**-1*g**-1'],name='cold',group=self)
        warm_filt = init_filter(cond_strs=['entropy_specific/>=/4.4e+8/erg*K**-1*g**-1','entropy_specific/<=/23.2e+8/erg*K**-1*g**-1'],name='warm',group=self)
        hot_filt  = init_filter(cond_strs=['entropy_specific/>/23.2e+8/erg*K**-1*g**-1'],name='hot',group=self)
        filt = [all_filt,cold_filt,warm_filt,hot_filt]

        # Since the number of global quanties that can be computed
        # from the gas data depends on the specific configuration
        # of a simulation, this needs to be determined at the start
        # See: ozy/sim_attributes.py/assign_attributes

        nvar_metals,nvar_magnetic,nvar_crs,nvar_rt = 0,0,0,0
        quantity_names = []
        do_binning = []
        weight_names = ['cumulative','density','mass','volume']
        if self.obj.simulation.physics['hydro']:
            quantity_names += ['mass','density','temperature',
                                'ang_momentum_x','ang_momentum_y',
                                'ang_momentum_z','thermal_energy',
                                'thermal_energy_specific',
                                'grav_therpfrsphere',
                                'sigma']
            do_binning += [False,True,True,
                            False,False,False,
                            False,True,True,True]
            if not self.obj.simulation.physics['magnetic'] and not \
                self.obj.simulation.physics['cr']:
                quantity_names += ['eff_FK2']
                do_binning += [True]
        nvar_metals = len(quantity_names)
        if self.obj.simulation.physics['metals']:
            quantity_names += ['metallicity']
            do_binning += [True]
        nvar_magnetic = len(quantity_names)
        if self.obj.simulation.physics['magnetic']:
            quantity_names += ['magnetic_energy','magnetic_energy_specific',
                                'magnetic_magnitude']
            do_binning += [False,True,True]
            if not self.obj.simulation.physics['cr']:
                quantity_names += ['eff_FKmag']
                do_binning += [True]
        nvar_crs = len(quantity_names)
        if self.obj.simulation.physics['cr']:
            quantity_names += ['cr_energy','cr_energy_specific',
                                'grav_crpfrsphere',
                                'eff_FKmag','eff_FKmagnocr']
            do_binning += [False,True,True,True,True]
        nvar_rt = len(quantity_names)
        if self.obj.simulation.physics['rt']:
            quantity_names += ['xHII','xHeII','xHeIII']
            do_binning += [False,False,False]
                
        # Initialise Fortran derived type with attributes
        # This object hold the following attributes:
        # - nvars: number of variables
        # - nwvars:  number of variables for weighting
        # - varnames: names of variables
        # - wvarnames:  names of variables for weighting
        # - data: organised in numpy array of shape (nvars,nwvars,4)
        #           each of those 4 values are (final, min, max, sum of weights)
        pdf_bins = 100
        glob_attrs = amr_integrator.amr_region_attrs()
        glob_attrs.nvars = len(quantity_names)
        glob_attrs.nwvars = len(weight_names)
        glob_attrs.nfilter = len(phase_names) + 1
        amr_integrator.allocate_amr_regions_attrs(glob_attrs)
        for i in range(0, len(quantity_names)):
            glob_attrs.varnames.T.view('S128')[i] = quantity_names[i].ljust(128)
            glob_attrs.result[i].nbins = pdf_bins
            glob_attrs.result[i].nfilter = len(filt)
            glob_attrs.result[i].nwvars = len(weight_names)
            glob_attrs.result[i].varname = quantity_names[i]
            mybins = get_code_bins(self.obj,'gas/'+quantity_names[i],pdf_bins)
            glob_attrs.result[i].scaletype = mybins[1]
            stats_utils.allocate_pdf(glob_attrs.result[i])
            glob_attrs.result[i].bins = mybins[0]
            glob_attrs.result[i].do_binning = do_binning[i]
            for j in range(0, len(weight_names)):
                glob_attrs.result[i].wvarnames.T.view('S128')[j] = weight_names[j].ljust(128)
        
        for i in range(0, glob_attrs.nfilter):
            glob_attrs.filters[i] = filt[i]

        # Begin integration
        use_neigh = True
        amr_integrator.integrate_region(output_path,selected_reg,use_neigh,glob_attrs)

        # Assign results to galaxy object
        if self.obj.simulation.physics['hydro']:
            self.mass['gas'] = self.obj.quantity(glob_attrs.result[0].total[0,0,0], 'code_mass')
            self.gas_density['gas'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[1],0),'code_density')
            self.temperature['gas'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[2],0),'code_temperature')
            self.angular_mom['gas'] = self.obj.array(np.array([glob_attrs.result[3].total[0,0,0],glob_attrs.result[4].total[0,0,0],glob_attrs.result[5].total[0,0,0]]),
                                                             'code_mass*code_length*code_velocity')
            self.energies['thermal_energy'] = self.obj.quantity(glob_attrs.result[6].total[0,0,0], 'code_mass * code_velocity**2')
            self.energies['thermal_energy_specific'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[7],0),'code_specific_energy')
            self.pressure_support['grav_therpfrsphere'] = pdf_handler_to_stats(self.obj,glob_attrs.result[8],0)
            self.velocity_dispersion['gas_turbulent'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[9],0),'code_velocity')
            if not self.obj.simulation.physics['magnetic'] and not self.obj.simulation.physics['cr']:
                self.sf_efficiency['eff_FK2'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[10],0),'dimensionless')
            if self.obj.simulation.physics['metals']:
                self.metallicity['gas'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[11],0),'dimensionless')
            for i in range(0, len(phase_names)):
                self.mass['gas_'+phase_names[i]] = self.obj.quantity(glob_attrs.result[0].total[i+1,0,0], 'code_mass')
                print('Mass in %s gas is %.5f'%(phase_names[i],self.mass['gas_'+phase_names[i]].to('Msun')))
                self.gas_density[phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[1],i+1),'code_density')
                self.temperature[phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[2],i+1),'code_temperature')
                self.angular_mom['gas_'+phase_names[i]] = self.obj.array(np.array([glob_attrs.result[3].total[i+1,0,0],
                                                                                    glob_attrs.result[4].total[i+1,0,0],
                                                                                    glob_attrs.result[5].total[i+1,0,0]]),
                                                                        'code_mass*code_length*code_velocity')
                self.energies['thermal_energy_'+phase_names[i]] = self.obj.quantity(glob_attrs.result[6].total[i+1,0,0], 'code_mass * code_velocity**2')
                self.energies['thermal_energy_specific_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[7],i+1),'code_specific_energy')
                self.pressure_support['grav_therpfrsphere_'+phase_names[i]] = pdf_handler_to_stats(self.obj,glob_attrs.result[8],i+1)
                self.velocity_dispersion['gas_turbulent_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[9],i+1),'code_velocity')
                if self.obj.simulation.physics['metals']:
                    self.metallicity['gas_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[10],i+1),'dimensionless')
        else:
            self.mass['gas'] = self.obj.quantity(0.0, 'code_mass')
        
        if self.obj.simulation.physics['magnetic']:
            print('Computing magnetic energies')
            self.energies['magnetic_energy'] = self.obj.quantity(glob_attrs.result[nvar_magnetic].total[0,0,0], 'code_mass * code_velocity**2')
            self.energies['magnetic_energy_specific'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[nvar_magnetic+1],0),'code_specific_energy')
            self.magnetism['magnetic_magnitude'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[nvar_magnetic+2],0),'code_magnetic')
            if not self.obj.simulation.physics['cr']:
                self.sf_efficiency['eff_FKmag'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[nvar_magnetic+3],0),'dimensionless')
            for i in range(0,len(phase_names)):
                self.energies['magnetic_energy_'+phase_names[i]] = self.obj.quantity(glob_attrs.result[nvar_magnetic].total[i+1,0,0], 'code_mass * code_velocity**2')
                self.energies['magnetic_energy_specific_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[nvar_magnetic+1],i+1),'code_specific_energy')
                self.magnetism['magnetic_magnitude_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[nvar_magnetic+2],i+1),'code_magnetic')
        
        if self.obj.simulation.physics['cr']:
            print('Computing CR energies')
            self.energies['cr_energy'] = self.obj.quantity(glob_attrs.result[nvar_crs].total[0,0,0], 'code_mass * code_velocity**2')
            self.energies['cr_energy_specific'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[nvar_crs+1],0),'code_specific_energy')
            self.pressure_support['grav_crpfrsphere'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[nvar_crs+2],0),'dimensionless')
            self.sf_efficiency['eff_FKmag'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[nvar_crs+3],0),'dimensionless')
            self.sf_efficiency['eff_FKmagnocr'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[nvar_crs+4],0),'dimensionless')
            for i in range(0, len(phase_names)):
                self.energies['cr_energy_'+phase_names[i]] = self.obj.quantity(glob_attrs.result[nvar_crs].total[i+1,0,0], 'code_mass * code_velocity**2')
                self.energies['cr_energy_specific_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[nvar_crs+1],i+1),'code_specific_energy')
                self.pressure_support['grav_crpfrsphere_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[nvar_crs+2],i+1), 'dimensionless')
        if self.obj.simulation.physics['rt']:
            #TODO: This is not really correct, need to update for the format of AMR integrations
            print('Computing ionisation fractions')
            self.radiation['xHII'] = glob_attrs.data[0,nvar_rt,1:,:-1]
            self.radiation['xHeII'] = glob_attrs.data[0,nvar_rt+1,1:,:-1]
            self.radiation['xHeIII'] = glob_attrs.data[0,nvar_rt+2,1:,:-1]
            for i in range(0, len(phase_names)):
                self.radiation['xHII_'+phase_names[i]] = glob_attrs.data[i+1,nvar_rt,1:,:-1]
                self.radiation['xHeII_'+phase_names[i]] = glob_attrs.data[i+1,nvar_rt+1,1:,:-1]
                self.radiation['xHeIII_'+phase_names[i]] = glob_attrs.data[i+1,nvar_rt+2,1:,:-1]

    def _calculate_halo_gas_quantities(self):
        from ozy.utils import get_code_bins, pdf_handler_to_stats
        """Compute gas quantities: Metallicity, Temperature..."""
        output_path = self.obj.simulation.fullpath

        # Initialise region
        selected_reg = init_region(self,'sphere',rmin=(0.2,'rvir'),rmax=(1.0,'rvir'))

        # We do not want any particular filter, just simple integration will do
        all_filt = filtering.filter()

        # Since the number of global quanties that can be computed
        # from the gas data depends on the specific configuration
        # of a simulation, this needs to be determined at the start
        # See: ozy/sim_attributes.py/assign_attributes

        nvar_metals,nvar_magnetic,nvar_crs,nvar_rt = 0,0,0,0
        quantity_names = []
        do_binning = []
        weight_names = ['cumulative','density','mass','volume']
        if self.obj.simulation.physics['hydro']:
            quantity_names += ['mass','density','temperature',
                                'ang_momentum_x','ang_momentum_y',
                                'ang_momentum_z','thermal_energy',
                                'thermal_energy_specific',
                                'grav_therpfrsphere',
                                'sigma']
            do_binning += [False,True,True,
                            False,False,False,
                            False,True,True,
                            True]
        nvar_metals = len(quantity_names)
        if self.obj.simulation.physics['metals']:
            quantity_names += ['metallicity']
            do_binning += [True]
        nvar_magnetic = len(quantity_names)
        if self.obj.simulation.physics['magnetic']:
            quantity_names += ['magnetic_energy','magnetic_energy_specific',
                                'magnetic_magnitude']
            do_binning += [False,True,True]
        nvar_crs = len(quantity_names)
        if self.obj.simulation.physics['cr']:
            quantity_names += ['cr_energy','cr_energy_specific',
                                'grav_crpfrsphere']
            do_binning += [False,True,True]
        nvar_rt = len(quantity_names)
        if self.obj.simulation.physics['rt']:
            quantity_names += ['xHII','xHeII','xHeIII']
            do_binning += [False,False,False]
                
        # Initialise Fortran derived type with attributes
        # This object hold the following attributes:
        # - nvars: number of variables
        # - nwvars:  number of variables for weighting
        # - varnames: names of variables
        # - wvarnames:  names of variables for weighting
        # - data: organised in numpy array of shape (nvars,nwvars,4)
        #           each of those 4 values are (final, min, max, sum of weights)
        pdf_bins = 100
        glob_attrs = amr_integrator.amr_region_attrs()
        glob_attrs.nvars = len(quantity_names)
        glob_attrs.nwvars = len(weight_names)
        glob_attrs.nfilter = 1
        amr_integrator.allocate_amr_regions_attrs(glob_attrs)
        for i in range(0, len(quantity_names)):
            glob_attrs.varnames.T.view('S128')[i] = quantity_names[i].ljust(128)
            glob_attrs.result[i].nbins = pdf_bins
            glob_attrs.result[i].nfilter = 1
            glob_attrs.result[i].nwvars = len(weight_names)
            glob_attrs.result[i].varname = quantity_names[i]
            mybins = get_code_bins(self.obj,'gas/'+quantity_names[i],pdf_bins)
            glob_attrs.result[i].scaletype = mybins[1]
            stats_utils.allocate_pdf(glob_attrs.result[i])
            glob_attrs.result[i].bins = mybins[0]
            glob_attrs.result[i].do_binning = do_binning[i]
            for j in range(0, len(weight_names)):
                glob_attrs.result[i].wvarnames.T.view('S128')[j] = weight_names[j].ljust(128)
        
        
        glob_attrs.filters[0] = all_filt

        # Begin integration
        use_neigh = True
        amr_integrator.integrate_region(output_path,selected_reg,use_neigh,glob_attrs)

        print('Integrating gas quantities for the halo region...')

        # Assign results to galaxy object
        if self.obj.simulation.physics['hydro']:
            self.mass['halo_gas'] = self.obj.quantity(glob_attrs.result[0].total[0,0,0], 'code_mass')
            self.gas_density['halo_gas'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[1],0),'code_density')
            self.temperature['halo_gas'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[2],0),'code_temperature')
            self.angular_mom['halo_gas'] = self.obj.array(np.array([glob_attrs.result[3].total[0,0,0],glob_attrs.result[4].total[0,0,0],glob_attrs.result[5].total[0,0,0]]),
                                                             'code_mass*code_length*code_velocity')
            self.energies['halo_thermal_energy'] = self.obj.quantity(glob_attrs.result[6].total[0,0,0], 'code_mass * code_velocity**2')
            self.energies['halo_thermal_energy_specific'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[7],0),'code_specific_energy')
            self.pressure_support['halo_grav_therpfrsphere'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[8],0),'dimensionless')
            self.velocity_dispersion['halo_gas_turbulent'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[9],0),'code_velocity')
            if self.obj.simulation.physics['metals']:
                self.metallicity['halo_gas'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[10],0),'dimensionless')
        else:
            self.mass['halo_gas'] = self.obj.quantity(np.zeros((3,7)), 'code_mass')
        
        if self.obj.simulation.physics['magnetic']:
            self.energies['halo_magnetic_energy'] = self.obj.quantity(glob_attrs.result[nvar_magnetic].total[0,0,0], 'code_mass * code_velocity**2')
            self.energies['halo_magnetic_energy_specific'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[nvar_magnetic+1],0),'code_specific_energy')
            self.magnetism['halo_magnetic_magnitude'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[nvar_magnetic+2],0),'code_magnetic')
        
        if self.obj.simulation.physics['cr']:
            self.energies['halo_cr_energy'] = self.obj.quantity(glob_attrs.result[nvar_crs].total[0,0,0], 'code_mass * code_velocity**2')
            self.energies['halo_cr_energy_specific'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[nvar_crs+1],0),'code_specific_energy')
            self.pressure_support['grav_crpfrsphere'] = self.obj.array(pdf_handler_to_stats(self.obj,glob_attrs.result[nvar_crs+2],0),'dimensionless')

        if self.obj.simulation.physics['rt']:
            #TODO: This is not really correct, need to update for the format of AMR integrations
            print('Computing ionisation fractions')
            self.radiation['halo_xHII'] = glob_attrs.data[0,nvar_rt,1:,:-1]
            self.radiation['halo_xHeII'] = glob_attrs.data[0,nvar_rt+1,1:,:-1]
            self.radiation['halo_xHeIII'] = glob_attrs.data[0,nvar_rt+2,1:,:-1]

    def _calculate_velocity_dispersions(self):
        """Calculate velocity dispersions for the various components."""
        # TODO
        return

    def _calculate_substructure_tidal(self,tidal_method):
        """Compute the tidal radii of substructure found for this object."""
        from ozy.utils import tidal_radius
        myhalo = self.obj.halos[self.parent_halo_index]
        if myhalo != None:
            subs = myhalo.substructure_list
            for s in subs:
                if s.npart >= 1000:
                    s.radius[tidal_method] = tidal_radius(myhalo,s,method=tidal_method)
        
    def _get_tdyn(self):
        """
        Computes the dynamical time-scale tdyn as
        the time required for a test particle to complete
        one full orbit at 0.2 Rvir.

        tdyn = 2pi*sqrt(R^3/(GM))
        where we assume M = Mgas+Mstars*Mdm is measured 
        within 0.2 Rvir
        """
        from unyt import G
        
        Mtot = self.mass['dm'] + self.mass['baryon']
        r = 0.2*self.obj.halos[self.parent_halo_index].virial_quantities['radius']
        tdyn = 2*np.pi*np.sqrt(r**3/(G*Mtot))
        return tdyn
class Halo(Group):
    """Halo class which has different levels of the halo hierarchy."""
    obj_type = 'halo'
    def __init__(self, obj):
        super(Halo, self).__init__(obj)
        self.spin = 0
        self.type               = 'halo'
        self.galaxies           = []
        self.central_galaxy     = None
        self.satellite_galaxies = []
        self.galaxy_index_list = []
        
    @property
    def substructure_list(self):
        subs = []
        if self.nextsub == 0:
            print('This halo does not seem to have a substructure assigned!')
            return subs
        haloIDs = [i.ID for i in self.obj.halos]
        nexti = self.nextsub
        while nexti != -1:
            if len(np.where(nexti == haloIDs)[0]) != 0:
                nextindex = np.where(nexti == haloIDs)[0][0]
                subs.append(self.obj.halos[nextindex])
                nexti = self.obj.halos[nextindex].nextsub
            else:
                break
        return subs
        
    def _compute_subs(self,tidal_method):
        from ozy.utils  import tidal_radius
        subs = self.substructure_list
        for s in subs:
            if s.npart >= 1000:
                tr = tidal_radius(self,s,method=tidal_method)
                s.radius[tidal_method] = tr

class Region(Group):
    """Region class - for regions of cosmological boxes or isolated test cases"""
    obj_type = 'region'
    def __init__(self, obj):
        """Initialise basic properties of the group"""
        self.obj = obj
        self.type = 'region'
        self.units = 'code'
        self.position = np.array([0.5, 0.5, 0.5]) #Needs to be thought through
        self.velocity = np.array([0., 0., 0.])
        self.angular_mom = {}
        self.virial_quantities = {}
        self.energies = {}

    def set_default_properties(self, obj):
        """Short-cut to set default properties, as required by ozymandias"""
        self.position = obj.array([0.5, 0.5, 0.5], 'code_length')
        self.velocity = obj.array([0., 0., 0.], 'km/s')
        self.virial_quantities['radius'] = obj.quantity(0.5, 'code_length') # Largest possible spherical radius
        self.angular_mom['total'] = obj.array([0., 0., 1.], 'Msun * km * km / s')


def create_new_group(obj, grouptype):
    """Simple function to create a new instance of a specified :class:`group.Group`.
    """
    if grouptype == 'halo':
        return Halo(obj)
    elif grouptype == 'galaxy':
        return Galaxy(obj)
    elif grouptype=='region':
        return Region(obj)
