common_variables = dict(
    x = {'unit':'code_length','neighbour':False},
    y = {'unit':'code_length','neighbour':False},
    z = {'unit':'code_length','neighbour':False},
    d_euclid = {'unit':'code_length','neighbour':False},
    r_sphere = {'unit':'code_length','neighbour':False},
    theta_sphere = {'unit':'radian','neighbour':False},
    phi_sphere = {'unit':'radian','neighbour':False},
    r_cyl = {'unit':'code_length','neighbour':False},
    phi_cyl = {'unit':'radian','neighbour':False},
    mass = {'unit':'code_mass','neighbour':False},
    density = {'unit':'code_density','neighbour':False},
    velocity_x = {'unit':'code_velocity','neighbour':False},
    velocity_y = {'unit':'code_velocity','neighbour':False},
    velocity_z = {'unit':'code_velocity','neighbour':False},
    v_sphere_r = {'unit':'code_velocity','neighbour':False},
    v_sphere_phi = {'unit':'code_velocity','neighbour':False},
    v_sphere_theta = {'unit':'code_velocity','neighbour':False},
    v_cyl_r = {'unit':'code_velocity','neighbour':False},
    v_cyl_z = {'unit':'code_velocity','neighbour':False},
    v_cyl_phi = {'unit':'code_velocity','neighbour':False},
    v_magnitude = {'unit':'code_velocity','neighbour':False},
    v_squared = {'unit':'code_velocity*code_velocity','neighbour':False},
    metallicity = {'unit':'code_metallicity','neighbour':False},
    kinetic_energy = {'unit':'code_energy','neighbour':False},
    momentum_x = {'unit':'code_mass*code_velocity','neighbour':False},
    momentum_y = {'unit':'code_mass*code_velocity','neighbour':False},
    momentum_z = {'unit':'code_mass*code_velocity','neighbour':False},
    momentum = {'unit':'code_mass*code_velocity','neighbour':False},
    momentum_sphere_r = {'unit':'code_mass*code_velocity','neighbour':False},
    absmomentum_sphere_r = {'unit':'code_mass*code_velocity','neighbour':False},
    momentum_cyl_z = {'unit':'code_mass*code_velocity','neighbour':False},
    ang_momentum_x = {'unit':'code_mass*code_length*code_velocity','neighbour':False},
    ang_momentum_y = {'unit':'code_mass*code_length*code_velocity','neighbour':False},
    ang_momentum_z = {'unit':'code_mass*code_length*code_velocity','neighbour':False},
    ang_momentum = {'unit':'code_mass*code_length*code_velocity','neighbour':False},
    counts = {'unit':'dimensionless','neighbour':False},
    cumulative = {'unit':'dimensionless','neighbour':False}
    
)
grid_variables = dict(
    temperature = {'unit':'code_temperature','neighbour':False},
    entropy_specific = {'unit':'code_specific_entropy','neighbour':False},
    volume = {'unit':'code_length*code_length*code_length','neighbour':False},
    sound_speed = {'unit':'code_velocity','neighbour':False},
    rms_speed = {'unit':'code_velocity','neighbour':False},
    alfven_speed = {'unit':'code_velocity','neighbour':False},
    massflow_rate = {'unit':'code_mass/code_time','neighbour':False},
    massflow_rate_sphere_r = {'unit':'code_mass/code_time','neighbour':False},
    massflux_rate_sphere_r = {'unit':'code_mass/code_time/code_length/code_length','neighbour':False},
    B_left_x = {'unit':'code_magnetic','neighbour':False},
    B_left_y = {'unit':'code_magnetic','neighbour':False},
    B_left_z = {'unit':'code_magnetic','neighbour':False},
    B_right_x = {'unit':'code_magnetic','neighbour':False},
    B_right_y = {'unit':'code_magnetic','neighbour':False},
    B_right_z = {'unit':'code_magnetic','neighbour':False},
    thermal_pressure = {'unit':'code_pressure','neighbour':False},
    thermal_energy = {'unit':'code_energy','neighbour':False},
    thermal_energy_specific = {'unit':'code_specific_energy','neighbour':False},
    net_cooling = {'unit':'code_energy_density/code_time','neighbour':False},
    grad_therprsphere = {'unit':'code_pressure/code_length','neighbour':True},
    grad_therpz = {'unit':'code_pressure/code_length','neighbour':True},
    total_coolingtime = {'unit':'code_time','neighbour':True},
    magnetic_magnitude = {'unit':'code_magnetic','neighbour':False},
    magnetic_energy = {'unit':'code_energy','neighbour':False},
    magnetic_energy_specific = {'unit':'code_specific_energy','neighbour':False},
    magnetic_energy_density = {'unit':'code_energy_density','neighbour':False},
    magnetic_pressure = {'unit':'code_pressure','neighbour':False},
    cr_energy = {'unit':'code_energy','neighbour':False},
    cr_energy_density = {'unit':'code_energy_density','neighbour':False},
    cr_pressure = {'unit':'code_pressure','neighbour':False},
    grad_crp = {'unit':'code_pressure/code_length','neighbour':True},
    grad_crprsphere = {'unit':'code_pressure/code_length','neighbour':True},
    grad_crpz = {'unit':'code_pressure/code_length','neighbour':True},
    gradscale_crprsphere = {'unit':'code_length','neighbour':True},
    gradscale_crp = {'unit':'code_length','neighbour':True},
    cr_energy_specific = {'unit':'code_specific_energy','neighbour':False},
    diffusion_speed = {'unit':'code_velocity','neighbour':True},
    alfvendiff_ratio = {'unit':'dimensionless','neighbour':True},
    stheatcooling_ratio = {'unit':'dimensionless','neighbour':True},
    cr_temperature_eff = {'unit':'code_temperature','neighbour':False},
    cr_GH08heat = {'unit':'code_energy_density/code_time','neighbour':False},
    streaming_heating = {'unit':'code_energy_density/code_time','neighbour':True},
    grav_crpf = {'unit':'dimensionless','neighbour':True},
    grav_crpfz = {'unit':'dimensionless','neighbour':True},
    grav_therpfz = {'unit':'dimensionless','neighbour':True},
    grav_crpfrsphere = {'unit':'dimensionless','neighbour':True},
    grav_crpfrspherepos = {'unit':'dimensionless','neighbour':True},
    grav_crpfrsphereneg = {'unit':'dimensionless','neighbour':True},
    grav_therpfrsphere = {'unit':'dimensionless','neighbour':True},
    grav_therpfrspherepos = {'unit':'dimensionless','neighbour':True},
    grav_therpfrsphereneg = {'unit':'dimensionless','neighbour':True},
    grav_totpfrsphere = {'unit':'dimensionless','neighbour':True},
    grav_totpfrspherepos = {'unit':'dimensionless','neighbour':True},
    grav_totpfrsphereneg = {'unit':'dimensionless','neighbour':True},
    grav_magpfrsphere = {'unit':'dimensionless','neighbour':True},
    grav_magpfrspherepos = {'unit':'dimensionless','neighbour':True},
    grav_magpfrsphereneg = {'unit':'dimensionless','neighbour':True},
    grav_gz = {'unit':'code_specific_energy/code_length','neighbour':False},
    grav_frsphere = {'unit':'code_pressure/code_length','neighbour':False},
    grav_fz = {'unit':'code_pressure/code_length','neighbour':False},
    xHII = {'unit':'dimensionless','neighbour':False},
    xHeII = {'unit':'dimensionless','neighbour':False},
    xHeIII = {'unit':'dimensionless','neighbour':False},
    dust_density = {'unit':'code_density','neighbour':False},
    DTM = {'unit':'dimensionless','neighbour':False},
    eff_FKmag = {'unit':'dimensionless','neighbour':True},
    eff_FKmagnocr = {'unit':'dimensionless','neighbour':True},
    eff_FK2 = {'unit':'dimensionless','neighbour':True},
    neighbour_accuracy = {'unit':'dimensionless','neighbour':True},
    sigma = {'unit':'code_velocity','neighbour':True}
)
particle_variables = dict(
    age = {'unit':'Gyr','neighbour':False},
    tform = {'unit':'code_time','neighbour':False},
    sfr = {'unit':'code_mass/code_time','neighbour':False},
    sfr_density = {'unit':'code_mass/Gyr/code_length**3','neighbour':False},
    sfr_surface = {'unit':'code_density*code_velocity','neighbour':False},
    sdensity = {'unit':'code_density*code_length','neighbour':False}
)

basic_conv = dict(
    code_length = 'kpc',
    code_mass = 'Msun',
    code_density = 'g*cm**-3',
    code_velocity = 'km*s**-1',
    code_energy = 'erg',
    code_specific_energy = 'erg*g**-1',
    code_specific_energy_code_length = 'cm*s**-2',
    code_energy_density = 'erg*cm**-3',
    code_energy_density_code_time = 'erg*s**-1*cm**-3',
    code_pressure = 'erg*cm**-3',
    code_pressure_code_length = 'erg*cm**-4',
    code_specific_entropy = 'erg*K**-1*g**-1',
    code_density_code_velocity_code_velocity = 'erg*cm**-3',
    code_density_code_velocity = 'Msun*yr**-1*kpc**-2',
    code_density_code_length = 'Msun*kpc**-2',
    code_mass_code_time = 'Msun*yr**-1',
    dimensionless = 'dimensionless',
    radian ='radian',
    code_magnetic = 'gauss',
    code_temperature = 'K',
    code_metallicity = 'dimensionless',
    code_time = 'yr',
    Gyr = 'Gyr'
)

def get_code_units(varname):

    if varname in common_variables:
        unit = common_variables[varname]['unit']
    elif varname in grid_variables:
        unit = grid_variables[varname]['unit']
    elif varname in particle_variables:
        unit = particle_variables[varname]['unit']
    else:
        raise KeyError('Variable not found, check: '+str(varname))
    return unit
    
def check_need_neighbours(varname):
    if varname in common_variables:
        need = common_variables[varname]['neighbour']
    elif varname in grid_variables:
        need = grid_variables[varname]['neighbour']
    elif varname in particle_variables:
        need = particle_variables[varname]['neighbour']
    else:
        raise KeyError('Variable not found, check: '+str(varname))
    return need