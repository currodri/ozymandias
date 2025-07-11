import numpy as np
from .group import Group
from .utils import init_region,get_code_bins,init_filter_hydro,\
                        init_filter_part, structure_regions

def integrate_hydro(obj, group=None, filter=None, xmin=(0,'code_length'), 
                    xmax=(1,'code_length'), ymin=(0,'code_length'),
                     ymax=(1,'code_length'), zmin=(0,'code_length'), 
                     zmax=(1,'code_length'),rmin=(0,'code_length'),
                     rmax=(0.5,'code_length'),region_type='sphere',
                     myaxis=np.array([0,0,1.]), bulk_velocity=([0,0,0],'code_velocity'),
                     variables=['mass'], weights=['cumulative'], do_binning=[False],
                     verbose=False, remove_subs=False, lmax=100, lmin=1, pdf_bins=100):
    """
    Integrate the hydro properties of a region in the simulation.

    Parameters
    ----------
    obj : object
        The simulation object.
    group : object, optional
        The group object to be used for integration. If None, a new group will be created.
    filter : object, optional
        The filter object to be used for integration. If None, a new filter will be created.
    xmin : tuple, optional
        The minimum x-coordinate of the region. Default is (0,'code_length').
    xmax : tuple, optional
        The maximum x-coordinate of the region. Default is (1,'code_length').
    ymin : tuple, optional
        The minimum y-coordinate of the region. Default is (0,'code_length').
    ymax : tuple, optional
        The maximum y-coordinate of the region. Default is (1,'code_length').
    zmin : tuple, optional
        The minimum z-coordinate of the region. Default is (0,'code_length').
    zmax : tuple, optional
        The maximum z-coordinate of the region. Default is (1,'code_length').
    rmin : tuple, optional
        The minimum radius of the region. Default is (0,'code_length').
    rmax : tuple, optional
        The maximum radius of the region. Default is (0.5,'code_length').
    region_type : str, optional
        The type of region to be integrated. Default is 'sphere'.
    myaxis : numpy array, optional
        The axis of rotation for the region. Default is np.array([0,0,1.]).
    bulk_velocity : tuple, optional
        The bulk velocity of the region. Default is ([0,0,0],'code_velocity').
    variables : list, optional
        The variables to be integrated. Default is ['mass'].
    weights : list, optional
        The weights to be used for integration. Default is ['cumulative'].
    do_binning : list, optional
        A list of booleans indicating whether to bin the variables. Default is [False].
    recompute : bool, optional
        If True, recompute the integration even if it has been done before. Default is False.
    verbose : bool, optional
        If True, print verbose output. Default is False.
    remove_subs : bool, optional
        If True, remove substructures from the integration. Default is False.
    lmax : int, optional
        The maximum level of refinement for the integration. Default is 100.
    lmin : int, optional
        The minimum level of refinement for the integration. Default is 1.

    Returns
    -------
    glob_attrs : object
        The global attributes of the integrator after the integration is performed.

    """
    from amr2 import io_ramses, amr_integrator, stats_utils
    from .utils import check_need_gravity, check_need_neighbours

    if verbose:
        io_ramses.activate_verbose()
    
    # 1. Initialise region
    if group == None:
        group = Group(obj)
        xcentre = 0.5*(obj.quantity(xmax[0],str(xmax[1])) + obj.quantity(xmin[0],str(xmin[1])))
        ycentre = 0.5*(obj.quantity(ymax[0],str(ymax[1])) + obj.quantity(ymin[0],str(ymin[1])))
        zcentre = 0.5*(obj.quantity(zmax[0],str(zmax[1])) + obj.quantity(zmin[0],str(zmin[1])))
        group.position = obj.array([xcentre.d,ycentre.d,zcentre.d],'code_length')
        group.angular_mom['total'] = myaxis
        group.velocity = obj.array(bulk_velocity[0],bulk_velocity[1])
        selected_reg = init_region(group, region_type, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                                    zmin=zmin, zmax=zmax, rmin=rmin, rmax=rmax)
        remove_subs = False # When group is None, we cannot remove substructures
    else:
        selected_reg,enclosing_sphere_p,enclosing_sphere_r = init_region(group, region_type, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                                    zmin=zmin, zmax=zmax, rmin=rmin, rmax=rmax, return_enclosing_sphere=True)
    
    # 2. If asked to remove substructures, obtain the list of substructures
    remove_all = False
    if remove_subs == 'all':
        remove_all = True
        remove_subs = True
    if remove_all:
        subs = structure_regions(group, add_substructure=False, add_neighbours=False,
                                    add_intersections=True,position=enclosing_sphere_p,
                                    radius=enclosing_sphere_r)
        nsubs = len(subs)
    elif remove_subs:
        subs = structure_regions(group, add_substructure=True, add_neighbours=False,
                                    tidal_method='BT87_simple')
        nsubs = len(subs)
    else:
        nsubs = 0 

    # 3. Get filter if not given
    if filter == None:
        filter = [init_filter_hydro('none','none',group)]
    elif not isinstance(filter, list):
        filter = [filter]

    # 4. Initialise the Fortran derived type for the integrator attributes
    glob_attrs = amr_integrator.amr_region_attrs()
    glob_attrs.nvars = len(variables)
    glob_attrs.nwvars = len(weights)
    glob_attrs.nfilter = len(filter)
    glob_attrs.nsubs = nsubs
    amr_integrator.allocate_amr_regions_attrs(glob_attrs)

    # 5. Now fill up the variable and weight names
    for i in range(0, len(variables)):
        glob_attrs.varnames.T.view('S128')[i] = variables[i].ljust(128)
    for i in range(0, len(weights)):
        glob_attrs.wvarnames.T.view('S128')[i] = weights[i].ljust(128)
    
    # 6. Now fill up the filters
    for i in range(0, len(filter)):
        glob_attrs.filters[i] = filter[i]

    # 7. Now fill up the substructure regions
    for i in range(0, nsubs):
        glob_attrs.subs[i] = subs[i]

    # 8. Setup of the PDF handler
    glob_attrs.result.nbins = max(pdf_bins,1)
    glob_attrs.result.nvars = len(variables)
    glob_attrs.result.nwvars = len(weights)
    glob_attrs.result.nfilter = len(filter)
    stats_utils.allocate_pdf(glob_attrs.result)

    # 9. Now fill up the variable and weight names with the correct bins
    for i in range(0, len(variables)):
        mybins = get_code_bins(obj, 'gas', variables[i], pdf_bins)
        glob_attrs.result.varname.T.view('S128')[i] = variables[i].ljust(128)
        glob_attrs.result.scaletype.T.view('S128')[i] = mybins[1].ljust(128)
        glob_attrs.result.bins[:,i] = mybins[0]
        glob_attrs.result.do_binning[i] = do_binning[i]
        glob_attrs.result.zero_index[i] = mybins[2]
        glob_attrs.result.linthresh[i] = mybins[3]
    for i in range(0, len(weights)):
        glob_attrs.result.wvarnames.T.view('S128')[i] = weights[i].ljust(128)

    # 10. Determine if the variables or weights need neighbours or gravity files
    use_neighbours = False
    use_gravity = False
    for i in range(0, len(variables)):
        if check_need_neighbours(variables[i],'gas'):
            use_neighbours = True
        if check_need_gravity(variables[i],'gas'):
            use_gravity = True
        if use_gravity and use_neighbours:
            break

    # 11. PERFORM THE INTEGRATION
    output_path = obj.simulation.fullpath
    if obj.use_vardict:
        amr_integrator.integrate_region(output_path, selected_reg,use_neighbours,
                                        use_gravity, glob_attrs, lmax, lmin,
                                        obj.vardict)
    else:
        amr_integrator.integrate_region(output_path, selected_reg,use_neighbours,
                                        use_gravity, glob_attrs, lmax, lmin)
    
    return glob_attrs

def integrate_part(obj, group=None, filter=None, xmin=(0,'code_length'),
                    xmax=(1,'code_length'), ymin=(0,'code_length'),
                    ymax=(1,'code_length'), zmin=(0,'code_length'), 
                    zmax=(1,'code_length'),rmin=(0,'code_length'),
                    rmax=(0.5,'code_length'),region_type='sphere',
                    myaxis=np.array([0,0,1.]), bulk_velocity=([0,0,0],'code_velocity'),
                    variables=['mass'], weights=['cumulative'], do_binning=[False],
                    verbose=False, remove_subs=False, pdf_bins=100):
    """
    Integrate the particle properties of a region in the simulation.

    Parameters
    ----------
    obj : object
        The simulation object.
    group : object, optional
        The group object to be used for integration. If None, a new group will be created.
    filter : object, optional
        The filter object to be used for integration. If None, a new filter will be created.
    xmin : tuple, optional
        The minimum x-coordinate of the region. Default is (0,'code_length').
    xmax : tuple, optional
        The maximum x-coordinate of the region. Default is (1,'code_length').
    ymin : tuple, optional
        The minimum y-coordinate of the region. Default is (0,'code_length').
    ymax : tuple, optional
        The maximum y-coordinate of the region. Default is (1,'code_length').
    zmin : tuple, optional
        The minimum z-coordinate of the region. Default is (0,'code_length').
    zmax : tuple, optional
        The maximum z-coordinate of the region. Default is (1,'code_length').
    rmin : tuple, optional
        The minimum radius of the region. Default is (0,'code_length').
    rmax : tuple, optional
        The maximum radius of the region. Default is (0.5,'code_length').
    region_type : str, optional
        The type of region to be integrated. Default is 'sphere'.
    myaxis : numpy array, optional
        The axis of rotation for the region. Default is np.array([0,0,1.]).
    bulk_velocity : tuple, optional
        The bulk velocity of the region. Default is ([0,0,0],'code_velocity').
    variables : list, optional
        The variables to be integrated. Default is ['mass'].
    weights : list, optional
        The weights to be used for integration. Default is ['cumulative'].
    do_binning : list, optional
        A list of booleans indicating whether to bin the variables. Default is [False].
    verbose : bool, optional
        If True, print verbose output. Default is False.
    remove_subs : bool, optional
        If True, remove substructures from the integration. Default is False.
    pdf_bins : int, optional
        The number of bins for the PDF. Default is 100.
    Returns
    -------
    glob_attrs : object
        The global attributes of the integrator after the integration is performed.

    """
    from part2 import io_ramses, part_integrator, stats_utils

    if verbose:
        io_ramses.activate_verbose()

    # 1. Initialise region
    if group == None:
        group = Group(obj)
        xcentre = 0.5*(obj.quantity(xmax[0],str(xmax[1])) + obj.quantity(xmin[0],str(xmin[1])))
        ycentre = 0.5*(obj.quantity(ymax[0],str(ymax[1])) + obj.quantity(ymin[0],str(ymin[1])))
        zcentre = 0.5*(obj.quantity(zmax[0],str(zmax[1])) + obj.quantity(zmin[0],str(zmin[1])))
        group.position = obj.array([xcentre.d,ycentre.d,zcentre.d],'code_length')
        group.angular_mom['total'] = myaxis
        group.velocity = obj.array(bulk_velocity[0],bulk_velocity[1])
        selected_reg = init_region(group, region_type, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                                    zmin=zmin, zmax=zmax, rmin=rmin, rmax=rmax)
        remove_subs = False # When group is None, we cannot remove substructures
    else:
        selected_reg,enclosing_sphere_p,enclosing_sphere_r = init_region(group, region_type, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                                    zmin=zmin, zmax=zmax, rmin=rmin, rmax=rmax, return_enclosing_sphere=True)
        
    # 2. If asked to remove substructures, obtain the list of substructures
    remove_all = False
    if remove_subs == 'all':
        remove_all = True
        remove_subs = True
    if remove_all:
        subs = structure_regions(group, add_substructure=False, add_neighbours=False,
                                    add_intersections=True,position=enclosing_sphere_p,
                                    radius=enclosing_sphere_r)
        nsubs = len(subs)
    elif remove_subs:
        subs = structure_regions(group, add_substructure=True, add_neighbours=False,
                                    tidal_method='BT87_simple')
        nsubs = len(subs)
    else:
        nsubs = 0
    
    # 3. Get filter if not given
    if filter == None:
        filter = [init_filter_part('none','none',group)]
    elif not isinstance(filter, list):
        filter = [filter]

    # 4. Initialise the Fortran derived type for the integrator attributes
    glob_attrs = part_integrator.part_region_attrs()
    glob_attrs.nvars = len(variables)
    glob_attrs.nwvars = len(weights)
    glob_attrs.nfilter = len(filter)
    glob_attrs.nsubs = nsubs
    part_integrator.allocate_part_regions_attrs(glob_attrs)

    # 5. Now fill up the variable and weight names
    for i in range(0, len(variables)):
        glob_attrs.varnames.T.view('S128')[i] = variables[i].ljust(128)
    for i in range(0, len(weights)):
        glob_attrs.wvarnames.T.view('S128')[i] = weights[i].ljust(128)

    # 6. Now fill up the filters
    for i in range(0, len(filter)):
        glob_attrs.filters[i] = filter[i]
    
    # 7. Now fill up the substructure regions
    for i in range(0, nsubs):
        glob_attrs.subs[i] = subs[i]
    
    # 8. Setup of the PDF handler
    glob_attrs.result.nbins = max(pdf_bins,1)
    glob_attrs.result.nvars = len(variables)
    glob_attrs.result.nwvars = len(weights)
    glob_attrs.result.nfilter = len(filter)
    stats_utils.allocate_pdf(glob_attrs.result)

    # 9. Now fill up the variable and weight names with the correct bins
    for i in range(0, len(variables)):
        mybins = get_code_bins(obj, 'part', variables[i], pdf_bins)
        glob_attrs.result.varname.T.view('S128')[i] = variables[i].ljust(128)
        glob_attrs.result.scaletype.T.view('S128')[i] = mybins[1].ljust(128)
        glob_attrs.result.bins[:,i] = mybins[0]
        glob_attrs.result.do_binning[i] = do_binning[i]
        glob_attrs.result.zero_index[i] = mybins[2]
        glob_attrs.result.linthresh[i] = mybins[3]
    for i in range(0, len(weights)):
        glob_attrs.result.wvarnames.T.view('S128')[i] = weights[i].ljust(128)

    # 10. PERFORM THE INTEGRATION
    output_path = obj.simulation.fullpath
    if obj.use_part_vardict:
        part_integrator.integrate_region(output_path, selected_reg, 
                                         glob_attrs, obj.part_vardict, 
                                         obj.part_vartypes)
    else:
        part_integrator.integrate_region(output_path, selected_reg, 
                                         glob_attrs)
        
    return glob_attrs