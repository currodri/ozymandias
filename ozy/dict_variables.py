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
    B_left_x = 'code_magnetic',
    B_left_y = 'code_magnetic',
    B_left_z = 'code_magnetic',
    B_right_x = 'code_magnetic',
    B_right_y = 'code_magnetic',
    B_right_z = 'code_magnetic',
    thermal_pressure = 'code_pressure',
    thermal_energy = 'code_energy',
    magnetic_energy = 'code_energy',
    cr_energy = 'code_energy',
)
particle_variables = dict(
    age = 'code_time',
    tform = 'code_time',
)

def get_code_units(varname):

    unit = 'code_length'
    if varname in common_variables:
        unit = common_variables[varname]
    elif varname in grid_variables:
        unit = grid_variables[varname]
    elif varname in particle_variables:
        unit = particle_variables[varname]
    return unit