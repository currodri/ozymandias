import numpy as np
from ozy.variables_settings import (
    geometrical_variables,
    raw_gas_variables,
    derived_gas_variables,
    gravity_variables,
    raw_part_variables,
    derived_part_variables,
    star_variables,
)

# Compatibility variable groups built from the new settings dictionaries.
common_variables = set(geometrical_variables.keys())
grid_variables = set(raw_gas_variables.keys()) | set(derived_gas_variables.keys()) | set(gravity_variables.keys())
particle_variables = set(raw_part_variables.keys()) | set(derived_part_variables.keys()) | set(star_variables.keys())

def unigrid_amr(obj, group=None, filter=None, lmax =0, n=[100,100,100], vars=['gas/density'], xmin=(0,'code_length'), xmax=(1,'code_length'), ymin=(0,'code_length'),
                 ymax=(1,'code_length'), zmin=(0,'code_length'), zmax=(1,'code_length'), angmom = np.array([0,0,1]), symlog=True):
    
    from .utils import init_region,init_filter
    from .group import create_new_group

    from amr2 import export_amr

    output_path = obj.simulation.fullpath

    # Initialise region
    if group == None:
        fake_obj = create_new_group('galaxy')
        xcentre = 0.5*(obj.quantity(xmax[0],str(xmax[1])) + obj.quantity(xmin[0],str(xmin[1])))
        ycentre = 0.5*(obj.quantity(ymax[0],str(ymax[1])) + obj.quantity(ymin[0],str(ymin[1])))
        zcentre = 0.5*(obj.quantity(zmax[0],str(zmax[1])) + obj.quantity(zmin[0],str(zmin[1])))
        fake_obj.position = obj.array([xcentre.d,ycentre.d,zcentre.d],'code_length')
        fake_obj.angular_mom['total'] = angmom
        selected_reg = init_region(fake_obj, 'basic_cube', xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                                    zmin=zmin, zmax=zmax)
    else:
        selected_reg = init_region(group, 'basic_cube', xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                                    zmin=zmin, zmax=zmax)
    # Get filter if not given
    if filter == None:
        filt = init_filter('none','none',group)
    
    # Get supported variables
    hydrovars = []
    for var in vars:
        var_type = var.split('/')[0]
        var_name = var.split('/')[1]
        if var_type == 'gas':
            if var_name in common_variables or var_name in grid_variables:
                hydrovars.append(var_name)
            else:
                raise KeyError('This gas variable is not supported. Please check!: %s', var)
        else:
            raise KeyError('This function only support gas variables and you are asking for %s variables', var_type)
    
    # Now initialise and allocate AMR chunk

    gridchunk = export_amr.chunk_handler()
    gridchunk.nx,gridchunk.ny,gridchunk.nz = n[0],n[1],n[2]
    gridchunk.filt = filt
    gridchunk.nvars = len(hydrovars)

    export_amr.allocate_chunk_handler(gridchunk)
    for i in range(0, len(hydrovars)):
        gridchunk.varnames.T.view('S128')[i] = hydrovars[i].ljust(128)

    # And finally, extract the unigrid
    print('Extracting %i variables from the AMR grid...'%len(hydrovars))
    export_amr.get_unigrid(output_path,selected_reg,lmax,symlog,gridchunk)
    grid = np.copy(gridchunk.data)

    return grid

def export2skirt(obj, group=None, filter=None, lmax =0, var='dust_density', xmin=(0,'code_length'), xmax=(1,'code_length'), ymin=(0,'code_length'),
                 ymax=(1,'code_length'), zmin=(0,'code_length'), zmax=(1,'code_length'),rmin=(0,'kpc'),rmax=(0.5,'code_length'), angmom = np.array([0,0,1]), 
                 symlog=True, h=(30,'pc'), smoothmethod='constant',sedmethod='bruzual&charlot',recompute=False):
    
    import os
    import ozy
    from .utils import init_region,init_filter
    from .group import create_new_group

    from part2 import export_part
    from amr2 import export_amr

    output_path = obj.simulation.fullpath

    # Initialise region
    if group == None:
        fake_obj = create_new_group('galaxy')
        xcentre = 0.5*(obj.quantity(xmax[0],str(xmax[1])) + obj.quantity(xmin[0],str(xmin[1])))
        ycentre = 0.5*(obj.quantity(ymax[0],str(ymax[1])) + obj.quantity(ymin[0],str(ymin[1])))
        zcentre = 0.5*(obj.quantity(zmax[0],str(zmax[1])) + obj.quantity(zmin[0],str(zmin[1])))
        fake_obj.position = obj.array([xcentre.d,ycentre.d,zcentre.d],'code_length')
        fake_obj.angular_mom['total'] = angmom
        selected_reg = init_region(fake_obj, 'basic_cube', xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                                    zmin=zmin, zmax=zmax)
    else:
        selected_reg = init_region(group, 'sphere', rmin=rmin, rmax=rmax)

    # Convert smoothing length to the required units (kpc)
    h = obj.quantity(h[0],h[1])

    # Get filter if not given
    if filter == None:
        filt = init_filter('none','none',group)
    
    # Create name for output files
    outid = output_path.split('/')[-1][-5:]
    if group == None:
        outname = os.getcwd()+'/snap_'+outid+'_noneobj'
    else:
        outname = os.getcwd()+'/snap_'+outid+'_gal'+str(group.ID)

    # Perform particle export
    partfile = outname + '_stars.txt'
    if os.path.exists(partfile) and recompute:
        export_part.part2skirt(output_path,selected_reg,filt,h.to('kpc').d,smoothmethod,sedmethod,partfile)
    elif not os.path.exists(partfile):
        export_part.part2skirt(output_path,selected_reg,filt,h.to('kpc').d,smoothmethod,sedmethod,partfile)

    # Perform hydro export
    gasfile = outname + '_gas.txt'
    if os.path.exists(gasfile) and recompute:
        export_amr.amr2skirt(output_path,selected_reg,filt,var,gasfile)
    elif not os.path.exists(gasfile):
        export_amr.amr2skirt(output_path,selected_reg,filt,var,gasfile)

def export2disperse(obj, group=None, filter=None, xmin=(0,'code_length'), xmax=(1,'code_length'), ymin=(0,'code_length'),
                    ymax=(1,'code_length'), zmin=(0,'code_length'), zmax=(1,'code_length'),rmin=(0,'kpc'),rmax=(0.5,'code_length'), 
                    angmom = np.array([0,0,1]), probability=0.01,recompute=False):

    import os
    import ozy
    from .utils import init_region,init_filter
    from .group import create_new_group

    from part2 import export_part

    output_path = obj.simulation.fullpath

    # Initialise region
    if group == None:
        fake_obj = create_new_group('galaxy')
        xcentre = 0.5*(obj.quantity(xmax[0],str(xmax[1])) + obj.quantity(xmin[0],str(xmin[1])))
        ycentre = 0.5*(obj.quantity(ymax[0],str(ymax[1])) + obj.quantity(ymin[0],str(ymin[1])))
        zcentre = 0.5*(obj.quantity(zmax[0],str(zmax[1])) + obj.quantity(zmin[0],str(zmin[1])))
        fake_obj.position = obj.array([xcentre.d,ycentre.d,zcentre.d],'code_length')
        fake_obj.angular_mom['total'] = angmom
        selected_reg = init_region(fake_obj, 'basic_cube', xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                                    zmin=zmin, zmax=zmax)
    else:
        selected_reg = init_region(group, 'sphere', rmin=rmin, rmax=rmax)

    # Get filter if not given
    if filter == None:
        filt = init_filter('none','none',group)
    
    # Create name for output files
    outid = output_path.split('/')[-1][-5:]
    if group == None:
        outname = os.getcwd()+'/snap_'+outid+'_noneobj'
    else:
        outname = os.getcwd()+'/snap_'+outid+'_gal'+str(group.ID)

    # Perform particle export
    partfile = outname + '_dm.txt'
    if os.path.exists(partfile) and recompute:
        export_part.part2disperse(output_path,selected_reg,filt,probability,partfile)
    elif not os.path.exists(partfile):
        export_part.part2disperse(output_path,selected_reg,filt,probability,partfile)
        
def basicexport2txt(obj, group=None, filter=None, lmax =0, var='dust_density', xmin=(0,'code_length'), xmax=(1,'code_length'), ymin=(0,'code_length'),
                    ymax=(1,'code_length'), zmin=(0,'code_length'), zmax=(1,'code_length'),rmin=(0,'kpc'),rmax=(0.5,'code_length'), angmom = np.array([0,0,1]), 
                    symlog=True, h=(30,'pc'), smoothmethod='constant',sedmethod='bruzual&charlot',recompute=False):
    
    import os
    import ozy
    from .utils import init_region,init_filter
    from .group import create_new_group

    from part2 import export_part
    from amr2 import export_amr

    output_path = obj.simulation.fullpath

    # Initialise region
    if group == None:
        fake_obj = create_new_group('galaxy')
        xcentre = 0.5*(obj.quantity(xmax[0],str(xmax[1])) + obj.quantity(xmin[0],str(xmin[1])))
        ycentre = 0.5*(obj.quantity(ymax[0],str(ymax[1])) + obj.quantity(ymin[0],str(ymin[1])))
        zcentre = 0.5*(obj.quantity(zmax[0],str(zmax[1])) + obj.quantity(zmin[0],str(zmin[1])))
        fake_obj.position = obj.array([xcentre.d,ycentre.d,zcentre.d],'code_length')
        fake_obj.angular_mom['total'] = angmom
        selected_reg = init_region(fake_obj, 'basic_cube', xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                                    zmin=zmin, zmax=zmax)
    else:
        selected_reg = init_region(group, 'sphere', rmin=rmin, rmax=rmax)

    # Convert smoothing length to the required units (kpc)
    h = obj.quantity(h[0],h[1])

    # Get filter if not given
    if filter == None:
        filt = init_filter('none','none',group)
    
    # Create name for output files
    outid = output_path.split('/')[-1][-5:]
    if group == None:
        outname = os.getcwd()+'/snap_'+outid+'_noneobj'
    else:
        outname = os.getcwd()+'/snap_'+outid+'_gal'+str(group.ID)

    # Perform particle export
    partfile = outname + '_stars.txt'
    if os.path.exists(partfile) and recompute:
        export_part.part2file(output_path,selected_reg,filt,h.to('kpc').d,smoothmethod,partfile)
    elif not os.path.exists(partfile):
        export_part.part2file(output_path,selected_reg,filt,h.to('kpc').d,smoothmethod,partfile)

def part2list(obj, group=None, filter=None, vars=['x', 'y', 'z', 'vx', 'vy', 'vz', 'mass'],
              xmin=(0,'code_length'), xmax=(1,'code_length'), ymin=(0,'code_length'), ymax=(1,'code_length'),
              zmin=(0,'code_length'), zmax=(1,'code_length'), rmin=(0,'kpc'), rmax=(0.5,'code_length'),
              angmom=np.array([0,0,1])):
    """Extract selected particle variables in a region and return them as an array.

    Returns
    -------
    data : np.ndarray
        Array of shape (nvars, nparticles), where rows follow the same order as ``vars``.
    """
    from .utils import init_region,init_filter_part
    from .group import create_new_group
    from part2 import export_part

    output_path = obj.simulation.fullpath

    # Initialise region
    if group == None:
        fake_obj = create_new_group('galaxy')
        xcentre = 0.5*(obj.quantity(xmax[0],str(xmax[1])) + obj.quantity(xmin[0],str(xmin[1])))
        ycentre = 0.5*(obj.quantity(ymax[0],str(ymax[1])) + obj.quantity(ymin[0],str(ymin[1])))
        zcentre = 0.5*(obj.quantity(zmax[0],str(zmax[1])) + obj.quantity(zmin[0],str(zmin[1])))
        fake_obj.position = obj.array([xcentre.d,ycentre.d,zcentre.d],'code_length')
        fake_obj.angular_mom['total'] = angmom
        selected_reg = init_region(fake_obj, 'basic_cube', xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                                   zmin=zmin, zmax=zmax)
    else:
        selected_reg = init_region(group, 'sphere', rmin=rmin, rmax=rmax)

    # Get filter if not given
    if filter == None:
        filt = init_filter_part('none','none',group)
    else:
        filt = filter

    # Keep variable handling consistent with integrate_part: pass requested
    # particle variable names directly to the Fortran tools.
    pvars = list(vars)

    nvars = len(pvars)

    # Use a derived-type handler so f90wrap can expose results reliably.
    chunk = export_part.part_chunk_handler()
    chunk.nvars = nvars
    chunk.filt = filt
    export_part.allocate_part_chunk_handler(chunk)

    for i in range(0, nvars):
        chunk.varnames.T.view('S128')[i] = pvars[i].ljust(128)

    if obj.use_part_vardict:
        export_part.part2list(
            output_path,
            selected_reg,
            chunk,
            obj.part_vardict,
            obj.part_vartypes,
        )
    else:
        export_part.part2list(output_path, selected_reg, chunk)

    return np.copy(chunk.data[:, :chunk.npart])

