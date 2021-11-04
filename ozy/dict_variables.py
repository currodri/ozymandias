common_variables = dict(
    d_euclid = 'code_length',
    r_sphere = 'code_length',
    theta_sphere = 'radian',
    phi_sphere = 'radian',
    r_cyl = 'code_length',
    phi_cyl = 'radian',
    mass = 'code_mass',
    density = 'code_density',
    velocity_x = 'code_velocity',
    velocity_y = 'code_velocity',
    velocity_z = 'code_velocity',
    v_sphere_r = 'code_velocity',
    v_sphere_phi = 'code_velocity',
    v_sphere_theta = 'code_velocity',
    v_cyl_r = 'code_velocity',
    v_cyl_z = 'code_velocity',
    v_cyl_phi = 'code_velocity',
    metallicity = 'code_metallicity',
    kinetic_energy = 'code_energy',
    momentum_x = 'code_mass*code_velocity',
    momentum_y = 'code_mass*code_velocity',
    momentum_z = 'code_mass*code_velocity',
    momentum = 'code_mass*code_velocity',
    momentum_sphere_r = 'code_mass*code_velocity',
    ang_momentum_x = 'code_mass*code_length*code_velocity',
    ang_momentum_y = 'code_mass*code_length*code_velocity',
    ang_momentum_z = 'code_mass*code_length*code_velocity',
    ang_momentum = 'code_mass*code_length*code_velocity',
    counts = 'dimensionless',
    cumulative = 'dimensionless'
    
)
grid_variables = dict(
    temperature = 'code_temperature',
    volume = 'code_volume',
    massflow_rate = 'code_mass/code_time',
    massflow_rate_sphere_r = 'code_mass/code_time',
    B_left_x = 'code_magnetic',
    B_left_y = 'code_magnetic',
    B_left_z = 'code_magnetic',
    B_right_x = 'code_magnetic',
    B_right_y = 'code_magnetic',
    B_right_z = 'code_magnetic',
    thermal_pressure = 'code_pressure',
    thermal_energy = 'code_energy',
    thermal_energy_specific = 'code_specific_energy',
    magnetic_magnitude = 'code_magnetic',
    magnetic_energy = 'code_energy',
    magnetic_energy_specific = 'code_specific_energy',
    magnetic_energy_density = 'code_density*code_velocity*code_velocity',
    magnetic_pressure = 'code_pressure',
    cr_energy = 'code_energy',
    cr_energy_density = 'code_density*code_velocity*code_velocity',
    cr_pressure = 'code_pressure',
    cr_energy_specific = 'code_specific_energy',
    xHII = 'dimensionless',
    xHeII = 'dimensionless',
    xHeIII = 'dimensionless'
)
particle_variables = dict(
    age = 'Gyr',
    tform = 'code_time',
    sfr = 'code_mass/code_time',
    sfr_density = 'code_mass/Gyr/code_length**3',
    sfr_surface = 'code_density*code_velocity',
    sdensity = 'code_density*code_length'
)

basic_conv = dict(
    code_length = 'kpc',
    code_mass = 'Msun',
    code_density = 'g*cm**-3',
    code_velocity = 'km*s**-1',
    code_energy = 'erg',
    code_specific_energy = 'erg*g**-1',
    code_energy_density = 'erg*cm**-3',
    code_density_code_velocity_code_velocity = 'erg*cm**-3',
    code_density_code_velocity = 'Msun*yr**-1*kpc**-2',
    code_density_code_length = 'Msun*kpc**-2',
    dimensionless = 'dimensionless',
    radian = 'radian',
    code_magnetic = 'gauss',
    code_temperature = 'K',
    code_metallicity = 'dimensionless',
    code_time = 'yr',
    Gyr = 'Gyr'
)

def get_code_units(varname):

    if varname in common_variables:
        unit = common_variables[varname]
    elif varname in grid_variables:
        unit = grid_variables[varname]
    elif varname in particle_variables:
        unit = particle_variables[varname]
    else:
        raise KeyError('Variable not found, check: '+str(varname))
    return unit