import numpy as np
from pprint import pprint
from .amr.amr2_pkg import amr_integrator,stats_utils
from .part.part2_pkg import part_integrator
from .variables_settings import MINIMUM_DM_PER_HALO, MINIMUM_STARS_PER_GALAXY
from .utils import init_filter_hydro,init_filter_part

grouptypes = dict(
    halo='halos',
    galaxy='galaxies',
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
        self._index   = -1
        self.level    = -1
        self.host     = -1
        self.hostsub  = -1
        self.nsub     = -1
        self.nextsub  = -1
        self.npart    = 0
        self.tstep    = -1
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
        if self.obj_type == 'halo' and self.npart < MINIMUM_DM_PER_HALO:
            return False
        elif self.obj_type == 'galaxy' and self.npart < MINIMUM_STARS_PER_GALAXY:
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
        self.nstar = 0
        self.parent_halo_index = -1
        self.gas_density = {}
        self.sfr     = {}
        self.energies = {}
        self.metallicity = {}
        self.temperature = {}
        self.radiation = {}
        self.magnetism = {}
        self.dust = {}
        self.sf_efficiency = {}
        self.pressure_support = {}
        self.velocity_dispersion = {}
        self.stellar_properties = {}
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

        empty_array = np.zeros((4,7)).astype(np.float64)

        self.mass['dm'] = self.obj.quantity(0.0, 'code_mass')
        self.angular_mom['dm'] = self.obj.array(np.array([0,0,0]),
                                                             'code_mass*code_length*code_velocity')
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
        self.metallicity['stellar'] = self.obj.array(empty_array, 'dimensionless') # Mass-weighted average!
        self.angular_mom['stellar'] = self.obj.array(np.array([0,0,0]),
                                                             'code_mass*code_length*code_velocity')
        

        phase_names = ['cold','warm','hot']

        if self.obj.simulation.physics['hydro']:
            self.mass['gas'] = self.obj.quantity(0.0, 'code_mass')
            self.gas_density['gas'] = self.obj.array(empty_array, 'code_density')
            self.temperature['gas'] = self.obj.array(empty_array, 'code_temperature')
            self.angular_mom['gas'] = self.obj.array(np.array([0.0,0.0,0.0]).astype(np.float64),
                                                             'code_mass*code_length*code_velocity')
            self.energies['thermal_energy'] = self.obj.quantity(0.0, 'code_mass * code_velocity**2')
            self.energies['thermal_energy_specific'] = self.obj.array(empty_array, 'code_specific_energy')
            self.pressure_support['grav_therpfrsphere'] = self.obj.array(empty_array, 'dimensionless')
            self.velocity_dispersion['gas_turbulent'] = self.obj.array(empty_array, 'dimensionless')
            if not self.obj.simulation.physics['magnetic'] and not self.obj.simulation.physics['cr']:
                self.sf_efficiency['eff_FK2'] = self.obj.array(empty_array, 'dimensionless')
            if self.obj.simulation.physics['metals']:
                self.metallicity['gas'] = self.obj.array(empty_array, 'dimensionless')
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
            self.radiation['xHII'] = self.obj.array(empty_array,'dimensionless')
            self.radiation['xHeII'] = self.obj.array(empty_array,'dimensionless')
            self.radiation['xHeIII'] = self.obj.array(empty_array,'dimensionless')
            for i in range(0, len(phase_names)):
                self.radiation['xHII_'+phase_names[i]] = self.obj.array(empty_array,'dimensionless')
                self.radiation['xHeII_'+phase_names[i]] = self.obj.array(empty_array,'dimensionless')
                self.radiation['xHeIII_'+phase_names[i]] = self.obj.array(empty_array,'dimensionless')

        # Assign results to galaxy object
        phase_names = ['hot','warm_ionised','warm_neutral','cold']
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
            for i in range(0, len(phase_names)):
                self.mass['halo_gas_'+phase_names[i]] = self.obj.quantity(0.0, 'code_mass')
                self.gas_density['halo_gas_'+phase_names[i]] = self.obj.array(empty_array,'code_density')
                self.temperature['halo_gas_'+phase_names[i]] = self.obj.array(empty_array,'code_temperature')
                self.angular_mom['halo_gas_'+phase_names[i]] = self.obj.array(np.array([0.0,0.0,0.0]),
                                                                'code_mass*code_length*code_velocity')
                self.energies['halo_thermal_energy_'+phase_names[i]] = self.obj.quantity(0.0, 'code_mass * code_velocity**2')
                self.energies['halo_thermal_energy_specific_'+phase_names[i]] = self.obj.array(empty_array,'code_specific_energy')
                self.pressure_support['halo_grav_therpfrsphere_'+phase_names[i]] = self.obj.array(empty_array,'dimensionless')
                self.velocity_dispersion['halo_gas_turbulent_'+phase_names[i]] = self.obj.array(empty_array,'code_velocity')
            
        else:
            self.mass['halo_gas'] = self.obj.quantity(0.0, 'code_mass')
        
        if self.obj.simulation.physics['magnetic']:
            self.energies['halo_magnetic_energy'] = self.obj.quantity(0.0, 'code_mass * code_velocity**2')
            self.energies['halo_magnetic_energy_specific'] = self.obj.array(empty_array, 'code_specific_energy')
            self.magnetism['halo_magnetic_magnitude'] = self.obj.array(empty_array, 'code_magnetic')
            for i in range(0, len(phase_names)):
                self.energies['halo_magnetic_energy_'+phase_names[i]] = self.obj.quantity(0.0, 'code_mass * code_velocity**2')
                self.energies['halo_magnetic_energy_specific_'+phase_names[i]] = self.obj.array(empty_array,'code_specific_energy')
                self.magnetism['halo_magnetic_magnitude_'+phase_names[i]] = self.obj.array(empty_array,'code_magnetic')
        
        if self.obj.simulation.physics['cr']:
            self.energies['halo_cr_energy'] = self.obj.quantity(0.0, 'code_mass * code_velocity**2')
            self.energies['halo_cr_energy_specific'] = self.obj.array(empty_array, 'code_specific_energy')
            self.pressure_support['grav_crpfrsphere'] = self.obj.array(empty_array, 'dimensionless')
            for i in range(0, len(phase_names)):
                self.energies['halo_cr_energy_'+phase_names[i]] = self.obj.quantity(0.0, 'code_mass * code_velocity**2')
                self.energies['halo_cr_energy_specific_'+phase_names[i]] = self.obj.array(empty_array,'code_specific_energy')
                self.pressure_support['grav_crpfrsphere_'+phase_names[i]] = self.obj.array(empty_array,'dimensionless')


        if self.obj.simulation.physics['rt']:
            self.radiation['halo_xHII'] = self.obj.array(empty_array,'dimensionless')
            self.radiation['halo_xHeII'] = self.obj.array(empty_array,'dimensionless')
            self.radiation['halo_xHeIII'] = self.obj.array(empty_array,'dimensionless')
            for i in range(0, len(phase_names)):
                self.radiation['halo_xHII_'+phase_names[i]] = self.obj.array(empty_array,'dimensionless')
                self.radiation['halo_xHeII_'+phase_names[i]] = self.obj.array(empty_array,'dimensionless')
                self.radiation['halo_xHeIII_'+phase_names[i]] = self.obj.array(empty_array,'dimensionless')


    def _calculate_stardm_quantities(self):
        """Calculate total stellar and DM masses as well as
            some useful global quantities.
        """
        from .integrators import integrate_part
        from .utils import pdf_handler_to_stats

        # Setup the filters for all particles, stars and DM
        filters = [init_filter_part(cond_strs=['none'],name='all',obj=self.obj),
                  init_filter_part(cond_strs=['age/>/0.0/Gyr'],name='stars',obj=self.obj),
                    init_filter_part(cond_strs=['age/==/0.0/Gyr'],name='dm',obj=self.obj)]
        
        # Setup the variables and weights
        variables = ['mass','age','sfr_10','sfr_100',
                     'ang_momentum_x','ang_momentum_y','ang_momentum_z']
        weights = ['cumulative','mass','age']
        do_binning = [False,True,False,False,False,False,False]
        nvar_metallicity = len(variables)
        if self.obj.simulation.physics['metallicity']:
            variables += ['metallicity']
            do_binning += [True]
        nvar_metals = len(variables)
        metal_vars = []
        if self.obj.simulation.physics['metals']:
            for el,present in self.obj.simulation.elements.items():
                if present:
                    variables += [el+'_fraction']
                    metal_vars += [el+'_fraction']
                    do_binning += [True]

        # Begin integration
        glob_attrs = integrate_part(self.obj,group=self,rmin=(0.0,'rvir'),
                                    rmax=(0.2,'rvir'),region_type='sphere',
                                    filter=filters,variables=variables,
                                    weights=weights,do_binning=do_binning,
                                    verbose=False)

        # Stellar details
        self.mass['stellar'] = self.obj.quantity(glob_attrs.result.total[0,1,0,0], 'code_mass')
        self.nstar = glob_attrs.result.nvalues[0,1]
        print('Number of stellar particles in galaxy %s: %d,%d'%(self.ID,self.nstar,self.npart))
        self.stellar_properties['age'] = self.obj.array(pdf_handler_to_stats(self.obj,'part',glob_attrs.result,1,1), 'Gyr')
        self.sfr['10Myr'] = self.obj.quantity(glob_attrs.result.total[2,1,0,0],'Msun/yr')
        self.sfr['100Myr'] = self.obj.quantity(glob_attrs.result.total[3,1,0,0],'Msun/yr')
        self.angular_mom['stellar'] = self.obj.array(np.array([glob_attrs.result.total[4,1,0,0],glob_attrs.result.total[5,1,0,0],glob_attrs.result.total[6,1,0,0]]),
                                                             'code_mass*code_length*code_velocity')
        if self.obj.simulation.physics['metallicity']:
            self.stellar_properties['metallicity'] = self.obj.array(pdf_handler_to_stats(self.obj,'part',glob_attrs.result,nvar_metallicity,1), 'dimensionless') # Mass-weighted average!
        if self.obj.simulation.physics['metals']:
            for i in range(0, len(metal_vars)):
                self.stellar_properties[metal_vars[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'part',glob_attrs.result,nvar_metals+i,1), 'dimensionless')
        
        # DM details
        self.mass['dm'] = self.obj.quantity(glob_attrs.result.total[0,2,0,0], 'code_mass')
        self.ndm = glob_attrs.result.nvalues[0,2]
        print('Number of DM particles in galaxy %s (%s): %d,%d'%(self.ID,self.obj.halos[self.parent_halo_index].ID,self.ndm,self.obj.halos[self.parent_halo_index].npart))
        self.angular_mom['dm'] = self.obj.array(np.array([glob_attrs.result.total[4,2,0,0],glob_attrs.result.total[5,2,0,0],glob_attrs.result.total[6,2,0,0]]),
                                                             'code_mass*code_length*code_velocity')


    def _calculate_galaxy_gas_quantities(self):
        from .integrators import integrate_hydro
        from .utils import pdf_handler_to_stats
        """Compute gas quantities: Metallicity, Temperature..."""

        # Define phase filters (for ISM, based on entropy)
        all_filt = init_filter_hydro('none','none',obj=self.obj)
        phase_names = ['cold','warm','hot']
        cold_filt = init_filter_hydro(cond_strs=['entropy_specific/</4.4e+8/erg*K**-1*g**-1'],name='cold',obj=self.obj)
        warm_filt = init_filter_hydro(cond_strs=['entropy_specific/>=/4.4e+8/erg*K**-1*g**-1','entropy_specific/<=/23.2e+8/erg*K**-1*g**-1'],name='warm',obj=self.obj)
        hot_filt  = init_filter_hydro(cond_strs=['entropy_specific/>/23.2e+8/erg*K**-1*g**-1'],name='hot',obj=self.obj)
        filt = [all_filt,cold_filt,warm_filt,hot_filt]

        # Since the number of global quantities that can be computed
        # from the gas data depends on the specific configuration
        # of a simulation, this needs to be determined at the start
        # See: ozy/sim_attributes.py/assign_attributes

        nvar_metallicity,nvar_metals,nvar_dust,nvar_magnetic,nvar_crs,nvar_rt = 0,0,0,0,0,0
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
        nvar_metallicity = len(quantity_names)
        if self.obj.simulation.physics['metallicity']:
            quantity_names += ['metallicity']
            do_binning += [True]
        nvar_metals = len(quantity_names)
        met_vars = []
        if self.obj.simulation.physics['metals']:
            for el,present in self.obj.simulation.elements.items():
                if present:
                    quantity_names += [el+'_fraction']
                    met_vars += [el+'_fraction']
                    do_binning += [True]
        nvar_dust = len(quantity_names)
        if self.obj.simulation.physics['dust']:
            quantity_names += ['CSmall_mass','CLarge_mass',
                               'SilSmall_mass','SilLarge_mass',
                               'CSmall_density','CLarge_density',
                               'SilSmall_density','SilLarge_density',
                               'CSmall_fraction','CLarge_fraction',
                               'SilSmall_fraction','SilLarge_fraction',
                               'DTM','DTG','STL','CSR']
            do_binning += [True,True,True,True,
                            True,True,True,True,
                            False,False,False,False,
                            True,True,True,True]
            if self.obj.simulation.physics['PAHs']:
                quantity_names += ['qPAH',
                                   'PAHSmall_mass','PAHLarge_mass',
                                   'PAHSmall_density','PAHLarge_density',
                                   'PAHSmall_fraction','PAHLarge_fraction'
                                   ]
                do_binning += [True,True,True,True,True,False,False]
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
            do_binning += [True,True,True]
                
        # Begin integration
        glob_attrs = integrate_hydro(self.obj,group=self,rmin=(0.0,'rvir'),
                                     rmax=(0.2,'rvir'),region_type='sphere',
                                     filter=filt,variables=quantity_names,
                                     weights=weight_names,do_binning=do_binning,
                                     verbose=False)

        # Assign results to galaxy object
        if self.obj.simulation.physics['hydro']:
            self.mass['gas'] = self.obj.quantity(glob_attrs.result.total[0,0,0,0], 'code_mass')
            self.gas_density['gas'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,1,0),'code_density')
            self.temperature['gas'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,2,0),'code_temperature')
            self.angular_mom['gas'] = self.obj.array(np.array([glob_attrs.result.total[3,0,0,0],
                                                               glob_attrs.result.total[4,0,0,0],
                                                               glob_attrs.result.total[5,0,0,0]]),
                                                             'code_mass*code_length*code_velocity')
            self.energies['thermal_energy'] = self.obj.quantity(glob_attrs.result.total[6,0,0,0], 'code_mass * code_velocity**2')
            self.energies['thermal_energy_specific'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,7,0),'code_specific_energy')
            self.pressure_support['grav_therpfrsphere'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,8,0),'dimensionless')
            self.velocity_dispersion['gas_turbulent'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,9,0),'code_velocity')
            if not self.obj.simulation.physics['magnetic'] and not self.obj.simulation.physics['cr']:
                self.sf_efficiency['eff_FK2'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,10,0),'dimensionless')
            if self.obj.simulation.physics['metallicity']:
                self.metallicity['gas'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_metallicity,0),'dimensionless')
            if self.obj.simulation.physics['metals']:
                for i in range(0, len(met_vars)):
                    self.metallicity[met_vars[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,
                                                                                          nvar_metals+i,0),'dimensionless')
            for i in range(0, len(phase_names)):
                self.mass['gas_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[0,i+1,0,0], 'code_mass')
                print('Mass in %s gas is %.5f'%(phase_names[i],self.mass['gas_'+phase_names[i]].to('Msun')))
                self.gas_density[phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,1,i+1),'code_density')
                self.temperature[phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,2,i+1),'code_temperature')
                self.angular_mom['gas_'+phase_names[i]] = self.obj.array(np.array([glob_attrs.result.total[3,0,0,0],
                                                                                    glob_attrs.result.total[4,i+1,0,0],
                                                                                    glob_attrs.result.total[5,i+1,0,0]]),
                                                                        'code_mass*code_length*code_velocity')
                self.energies['thermal_energy_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[6,i+1,0,0], 'code_mass * code_velocity**2')
                self.energies['thermal_energy_specific_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,7,i+1),'code_specific_energy')
                self.pressure_support['grav_therpfrsphere_'+phase_names[i]] = pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,8,i+1)
                self.velocity_dispersion['gas_turbulent_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,9,i+1),'code_velocity')
                if self.obj.simulation.physics['metallicity']:
                    self.metallicity['gas_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_metallicity,i+1),'dimensionless')
                if self.obj.simulation.physics['metals']:
                    for j in range(0, len(met_vars)):
                        self.metallicity[met_vars[j]+'_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,
                                                                                                  nvar_metals+j,i+1),'dimensionless')
        else:
            self.mass['gas'] = self.obj.quantity(0.0, 'code_mass')

        if self.obj.simulation.physics['dust']:
            self.dust['CSmall_mass'] = self.obj.quantity(glob_attrs.result.total[nvar_dust,0,0,0], 'code_mass')
            self.dust['CLarge_mass'] = self.obj.quantity(glob_attrs.result.total[nvar_dust+1,0,0,0], 'code_mass')
            self.dust['SilSmall_mass'] = self.obj.quantity(glob_attrs.result.total[nvar_dust+2,0,0,0], 'code_mass')
            self.dust['SilLarge_mass'] = self.obj.quantity(glob_attrs.result.total[nvar_dust+3,0,0,0], 'code_mass')
            self.dust['CSmall_density'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+4,0),'code_density')
            self.dust['CLarge_density'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+5,0),'code_density')
            self.dust['SilSmall_density'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+6,0),'code_density')
            self.dust['SilLarge_density'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+7,0),'code_density')
            self.dust['CSmall_fraction'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+8,0),'dimensionless')
            self.dust['CLarge_fraction'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+9,0),'dimensionless')
            self.dust['SilSmall_fraction'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+10,0),'dimensionless')
            self.dust['SilLarge_fraction'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+11,0),'dimensionless')
            self.dust['DTM'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+12,0),'dimensionless')
            self.dust['DTG'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+13,0),'dimensionless')
            self.dust['STL'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+14,0),'dimensionless')
            self.dust['CSR'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+15,0),'dimensionless')
            if self.obj.simulation.physics['PAHs']:
                self.dust['qPAH'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+16,0),'dimensionless')
                self.dust['PAHSmall_mass'] = self.obj.quantity(glob_attrs.result.total[nvar_dust+17,0,0,0], 'code_mass')
                self.dust['PAHLarge_mass'] = self.obj.quantity(glob_attrs.result.total[nvar_dust+18,0,0,0], 'code_mass')
                self.dust['PAHSmall_density'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+19,0),'code_density')
                self.dust['PAHLarge_density'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+20,0),'code_density')
                self.dust['PAHSmall_fraction'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+21,0),'dimensionless')
                self.dust['PAHLarge_fraction'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+22,0),'dimensionless')
            for i in range(0, len(phase_names)):
                self.dust['CSmall_mass_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[nvar_dust,i+1,0,0], 'code_mass')
                self.dust['CLarge_mass_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[nvar_dust+1,i+1,0,0], 'code_mass')
                self.dust['SilSmall_mass_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[nvar_dust+2,i+1,0,0], 'code_mass')
                self.dust['SilLarge_mass_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[nvar_dust+3,i+1,0,0], 'code_mass')
                self.dust['CSmall_density_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+4,i+1),'code_density')
                self.dust['CLarge_density_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+5,i+1),'code_density')
                self.dust['SilSmall_density_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+6,i+1),'code_density')
                self.dust['SilLarge_density_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+7,i+1),'code_density')
                self.dust['CSmall_fraction_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+8,i+1),'dimensionless')
                self.dust['CLarge_fraction_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+9,i+1),'dimensionless')
                self.dust['SilSmall_fraction_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+10,i+1),'dimensionless')
                self.dust['SilLarge_fraction_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+11,i+1),'dimensionless')
                self.dust['DTM_'+phase_names[i]] = pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+12,i+1)
                self.dust['DTG_'+phase_names[i]] = pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+13,i+1)
                self.dust['STL_'+phase_names[i]] = pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+14,i+1)
                self.dust['CSR_'+phase_names[i]] = pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+15,i+1)
                if self.obj.simulation.physics['PAHs']:
                    self.dust['qPAH_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+16,i+1),'dimensionless')
                    self.dust['PAHSmall_mass_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[nvar_dust+17,i+1,0,0], 'code_mass')
                    self.dust['PAHLarge_mass_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[nvar_dust+18,i+1,0,0], 'code_mass')
                    self.dust['PAHSmall_density_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+19,i+1),'code_density')
                    self.dust['PAHLarge_density_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+20,i+1),'code_density')
                    self.dust['PAHSmall_fraction_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+21,i+1),'dimensionless')
                    self.dust['PAHLarge_fraction_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+22,i+1),'dimensionless')
        
        if self.obj.simulation.physics['magnetic']:
            self.energies['magnetic_energy'] = self.obj.quantity(glob_attrs.result.total[nvar_magnetic,0,0,0], 'code_mass * code_velocity**2')
            self.energies['magnetic_energy_specific'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_magnetic+1,0),'code_specific_energy')
            self.magnetism['magnetic_magnitude'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_magnetic+2,0),'code_magnetic')
            if not self.obj.simulation.physics['cr']:
                self.sf_efficiency['eff_FKmag'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_magnetic+3,0),'dimensionless')
            for i in range(0,len(phase_names)):
                self.energies['magnetic_energy_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[nvar_magnetic,i+1,0,0], 'code_mass * code_velocity**2')
                self.energies['magnetic_energy_specific_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_magnetic+1,i+1),'code_specific_energy')
                self.magnetism['magnetic_magnitude_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_magnetic+2,i+1),'code_magnetic')
        
        if self.obj.simulation.physics['cr']:
            self.energies['cr_energy'] = self.obj.quantity(glob_attrs.result.total[nvar_crs,0,0,0], 'code_mass * code_velocity**2')
            self.energies['cr_energy_specific'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_crs+1,0),'code_specific_energy')
            self.pressure_support['grav_crpfrsphere'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_crs+2,0),'dimensionless')
            self.sf_efficiency['eff_FKmag'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_crs+3,0),'dimensionless')
            self.sf_efficiency['eff_FKmagnocr'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_crs+4,0),'dimensionless')
            for i in range(0, len(phase_names)):
                self.energies['cr_energy_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[nvar_crs,i+1,0,0], 'code_mass * code_velocity**2')
                self.energies['cr_energy_specific_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_crs+1,i+1),'code_specific_energy')
                self.pressure_support['grav_crpfrsphere_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_crs+2,i+1), 'dimensionless')
        if self.obj.simulation.physics['rt']:
            self.radiation['xHII'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_rt,0),'dimensionless')
            self.radiation['xHeII'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_rt+1,0),'dimensionless')
            self.radiation['xHeIII'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_rt+2,0),'dimensionless')
            for i in range(0, len(phase_names)):
                self.radiation['xHII_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_rt,i+1),'dimensionless')
                self.radiation['xHeII_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_rt+1,i+1),'dimensionless')
                self.radiation['xHeIII_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_rt+2,i+1),'dimensionless')

    def _calculate_halo_gas_quantities(self):
        from .integrators import integrate_hydro
        from .utils import pdf_handler_to_stats
        """Compute gas quantities: Metallicity, Temperature..."""

        # Define phase filters (for CGM)
        all_filt = init_filter_hydro('none','none',obj=self.obj)
        phase_names = ['hot','warm_ionised','warm_neutral','cold']
        hot = init_filter_hydro(cond_strs=['temperature/>/1e5/K'],name='hot',obj=self.obj)
        warm_ionised = init_filter_hydro(cond_strs=['temperature/</1e5/K','temperature/>/9e3/K'],name='warm_ionised',obj=self.obj)
        warm_neutral = init_filter_hydro(cond_strs=['temperature/</9e3/K','temperature/>/1e3/K'],name='warm_neutral',obj=self.obj)
        cold = init_filter_hydro(cond_strs=['temperature/</1e3/K'],name='cold',obj=self.obj)
        filt = [all_filt,hot,warm_ionised,warm_neutral,cold]

        # Since the number of global quanties that can be computed
        # from the gas data depends on the specific configuration
        # of a simulation, this needs to be determined at the start
        # See: ozy/sim_attributes.py/assign_attributes

        nvar_metallicity,nvar_metals,nvar_dust,nvar_magnetic,nvar_crs,nvar_rt = 0,0,0,0,0,0
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
        nvar_metallicity = len(quantity_names)
        if self.obj.simulation.physics['metals']:
            quantity_names += ['metallicity']
            do_binning += [True]
        nvar_metals = len(quantity_names)
        met_vars = []
        if self.obj.simulation.physics['metals']:
            for el,present in self.obj.simulation.elements.items():
                if present:
                    quantity_names += [el+'_fraction']
                    met_vars += [el+'_fraction']
                    do_binning += [True]
        nvar_dust = len(quantity_names)
        if self.obj.simulation.physics['dust']:
            quantity_names += ['CSmall_mass','CLarge_mass',
                               'SilSmall_mass','SilLarge_mass',
                               'CSmall_density','CLarge_density',
                               'SilSmall_density','SilLarge_density',
                               'CSmall_fraction','CLarge_fraction',
                               'SilSmall_fraction','SilLarge_fraction',
                               'DTM','DTG','STL','CSR']
            do_binning += [True,True,True,True,
                            True,True,True,True,
                            False,False,False,False,
                            True,True,True,True]
            if self.obj.simulation.physics['PAHs']:
                quantity_names += ['qPAH',
                                   'PAHSmall_mass','PAHLarge_mass',
                                   'PAHSmall_density','PAHLarge_density',
                                   'PAHSmall_fraction','PAHLarge_fraction'
                                   ]
                do_binning += [True,True,True,True,True,False,False]
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
            do_binning += [True,True,True]
                
        # Begin integration
        glob_attrs = integrate_hydro(self.obj,group=self,rmin=(0.2,'rvir'),
                                     rmax=(1.0,'rvir'),region_type='sphere',
                                     filter=filt,variables=quantity_names,
                                     weights=weight_names,do_binning=do_binning)

        print('Integrating gas quantities for the halo region...')

        # Assign results to galaxy object
        if self.obj.simulation.physics['hydro']:
            self.mass['halo_gas'] = self.obj.quantity(glob_attrs.result.total[0,0,0,0], 'code_mass')
            self.gas_density['halo_gas'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,1,0),'code_density')
            self.temperature['halo_gas'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,2,0),'code_temperature')
            self.angular_mom['halo_gas'] = self.obj.array(np.array([glob_attrs.result.total[3,0,0,0],
                                                                    glob_attrs.result.total[4,0,0,0],
                                                                    glob_attrs.result.total[5,0,0,0]]),
                                                             'code_mass*code_length*code_velocity')
            self.energies['halo_thermal_energy'] = self.obj.quantity(glob_attrs.result.total[6,0,0,0], 'code_mass * code_velocity**2')
            self.energies['halo_thermal_energy_specific'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,7,0),'code_specific_energy')
            self.pressure_support['halo_grav_therpfrsphere'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,8,0),'dimensionless')
            self.velocity_dispersion['halo_gas_turbulent'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,9,0),'code_velocity')
            if self.obj.simulation.physics['metallicity']:
                self.metallicity['halo_gas'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_metallicity,0),'dimensionless')
            if self.obj.simulation.physics['metals']:
                for i in range(0, len(met_vars)):
                    self.metallicity[met_vars[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,
                                                                                          nvar_metals+i,0),'dimensionless')
            for i in range(0, len(phase_names)):
                self.mass['halo_gas_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[0,i+1,0,0], 'code_mass')
                self.gas_density['halo_gas_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,1,i+1),'code_density')
                self.temperature['halo_gas_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,2,i+1),'code_temperature')
                self.angular_mom['halo_gas_'+phase_names[i]] = self.obj.array(np.array([glob_attrs.result.total[3,i+1,0,0],
                                                                                        glob_attrs.result.total[4,i+1,0,0],
                                                                                        glob_attrs.result.total[5,i+1,0,0]]),
                                                                'code_mass*code_length*code_velocity')
                self.energies['halo_thermal_energy_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[6,i+1,0,0], 'code_mass * code_velocity**2')
                self.energies['halo_thermal_energy_specific_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,7,i+1),'code_specific_energy')
                self.pressure_support['halo_grav_therpfrsphere_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,8,i+1),'dimensionless')
                self.velocity_dispersion['halo_gas_turbulent_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,9,i+1),'code_velocity')
                if self.obj.simulation.physics['metallicity']:
                    self.metallicity['halo_gas_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_metallicity,i+1),'dimensionless')
                if self.obj.simulation.physics['metals']:
                    for j in range(0, len(met_vars)):
                        self.metallicity[met_vars[j]+'_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,
                                                                                                  nvar_metals+j,i+1),'dimensionless')
        else:
            self.mass['halo_gas'] = self.obj.quantity(np.zeros((3,7)), 'code_mass')
        
        if self.obj.simulation.physics['dust']:
            self.dust['halo_CSmall_mass'] = self.obj.quantity(glob_attrs.result.total[nvar_dust,0,0,0], 'code_mass')
            self.dust['halo_CLarge_mass'] = self.obj.quantity(glob_attrs.result.total[nvar_dust+1,0,0,0], 'code_mass')
            self.dust['halo_SilSmall_mass'] = self.obj.quantity(glob_attrs.result.total[nvar_dust+2,0,0,0], 'code_mass')
            self.dust['halo_SilLarge_mass'] = self.obj.quantity(glob_attrs.result.total[nvar_dust+3,0,0,0], 'code_mass')
            self.dust['halo_CSmall_density'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+4,0),'code_density')
            self.dust['halo_CLarge_density'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+5,0),'code_density')
            self.dust['halo_SilSmall_density'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+6,0),'code_density')
            self.dust['halo_SilLarge_density'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+7,0),'code_density')
            self.dust['halo_CSmall_fraction'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+8,0),'dimensionless')
            self.dust['halo_CLarge_fraction'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+9,0),'dimensionless')
            self.dust['halo_SilSmall_fraction'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+10,0),'dimensionless')
            self.dust['halo_SilLarge_fraction'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+11,0),'dimensionless')
            self.dust['halo_DTM'] = pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+12,0)
            self.dust['halo_DTG'] = pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+13,0)
            self.dust['halo_STL'] = pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+14,0)
            self.dust['halo_CSR'] = pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+15,0)
            if self.obj.simulation.physics['PAHs']:
                self.dust['halo_qPAH'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+16,0),'dimensionless')
                self.dust['halo_PAHSmall_mass'] = self.obj.quantity(glob_attrs.result.total[nvar_dust+17,0,0,0], 'code_mass')
                self.dust['halo_PAHLarge_mass'] = self.obj.quantity(glob_attrs.result.total[nvar_dust+18,0,0,0], 'code_mass')
                self.dust['halo_PAHSmall_density'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+19,0),'code_density')
                self.dust['halo_PAHLarge_density'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+20,0),'code_density')
                self.dust['halo_PAHSmall_fraction'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+21,0),'dimensionless')
                self.dust['halo_PAHLarge_fraction'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+22,0),'dimensionless')
            for i in range(0, len(phase_names)):
                self.dust['halo_CSmall_mass_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[nvar_dust,i+1,0,0], 'code_mass')
                self.dust['halo_CLarge_mass_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[nvar_dust+1,i+1,0,0], 'code_mass')
                self.dust['halo_SilSmall_mass_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[nvar_dust+2,i+1,0,0], 'code_mass')
                self.dust['halo_SilLarge_mass_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[nvar_dust+3,i+1,0,0], 'code_mass')
                self.dust['halo_CSmall_density_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+4,i+1),'code_density')
                self.dust['halo_CLarge_density_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+5,i+1),'code_density')
                self.dust['halo_SilSmall_density_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+6,i+1),'code_density')
                self.dust['halo_SilLarge_density_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+7,i+1),'code_density')
                self.dust['halo_CSmall_fraction_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+8,i+1),'dimensionless')
                self.dust['halo_CLarge_fraction_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+9,i+1),'dimensionless')
                self.dust['halo_SilSmall_fraction_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+10,i+1),'dimensionless')
                self.dust['halo_SilLarge_fraction_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+11,i+1),'dimensionless')
                self.dust['halo_DTM_'+phase_names[i]] = pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+12,i+1)
                self.dust['halo_DTG_'+phase_names[i]] = pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+13,i+1)
                self.dust['halo_STL_'+phase_names[i]] = pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+14,i+1)
                self.dust['halo_CSR_'+phase_names[i]] = pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+15,i+1)
                if self.obj.simulation.physics['PAHs']:
                    self.dust['halo_qPAH_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+16,i+1),'dimensionless')
                    self.dust['halo_PAHSmall_mass_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[nvar_dust+17,i+1,0,0], 'code_mass')
                    self.dust['halo_PAHLarge_mass_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[nvar_dust+18,i+1,0,0], 'code_mass')
                    self.dust['halo_PAHSmall_density_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+19,i+1),'code_density')
                    self.dust['halo_PAHLarge_density_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+20,i+1),'code_density')
                    self.dust['halo_PAHSmall_fraction_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+21,i+1),'dimensionless')
                    self.dust['halo_PAHLarge_fraction_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_dust+22,i+1),'dimensionless')

        if self.obj.simulation.physics['magnetic']:
            self.energies['halo_magnetic_energy'] = self.obj.quantity(glob_attrs.result.total[nvar_magnetic,0,0,0], 'code_mass * code_velocity**2')
            self.energies['halo_magnetic_energy_specific'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_magnetic+1,0),'code_specific_energy')
            self.magnetism['halo_magnetic_magnitude'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_magnetic+2,0),'code_magnetic')
            for i in range(0, len(phase_names)):
                self.energies['halo_magnetic_energy_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[nvar_magnetic,i+1,0,0], 'code_mass * code_velocity**2')
                self.energies['halo_magnetic_energy_specific_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_magnetic+1,i+1),'code_specific_energy')
                self.magnetism['halo_magnetic_magnitude_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_magnetic+2,i+1),'code_magnetic')
                
        if self.obj.simulation.physics['cr']:
            self.energies['halo_cr_energy'] = self.obj.quantity(glob_attrs.result.total[nvar_crs,0,0,0], 'code_mass * code_velocity**2')
            self.energies['halo_cr_energy_specific'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_crs+1,0),'code_specific_energy')
            self.pressure_support['grav_crpfrsphere'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_crs+2,0),'dimensionless')
            for i in range(0, len(phase_names)):
                self.energies['halo_cr_energy_'+phase_names[i]] = self.obj.quantity(glob_attrs.result.total[nvar_crs,i+1,0,0], 'code_mass * code_velocity**2')
                self.energies['halo_cr_energy_specific_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_crs+1,i+1),'code_specific_energy')
                self.pressure_support['grav_crpfrsphere_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_crs+2,i+1),'dimensionless')

        if self.obj.simulation.physics['rt']:
            self.radiation['halo_xHII'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_rt,0),'dimensionless')
            self.radiation['halo_xHeII'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_rt+1,0),'dimensionless')
            self.radiation['halo_xHeIII'] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_rt+2,0),'dimensionless')
            for i in range(0, len(phase_names)):
                self.radiation['halo_xHII_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_rt,i+1),'dimensionless')
                self.radiation['halo_xHeII_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_rt+1,i+1),'dimensionless')
                self.radiation['halo_xHeIII_'+phase_names[i]] = self.obj.array(pdf_handler_to_stats(self.obj,'gas',glob_attrs.result,nvar_rt+2,i+1),'dimensionless')


    def _calculate_velocity_dispersions(self):
        """Calculate velocity dispersions for the various components."""
        # TODO
        return

    def _calculate_substructure_tidal(self,tidal_method):
        """Compute the tidal radii of substructure found for this object."""
        from .utils import tidal_radius
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
        self.ndm               = 0
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
        from .utils  import tidal_radius
        subs = self.substructure_list
        for s in subs:
            if s.npart >= 1000:
                tr = tidal_radius(self,s,method=tidal_method)
                s.radius[tidal_method] = tr

def create_new_group(obj, grouptype):
    """Simple function to create a new instance of a specified :class:`group.Group`.
    """
    if grouptype == 'halo':
        return Halo(obj)
    elif grouptype == 'galaxy':
        return Galaxy(obj)
