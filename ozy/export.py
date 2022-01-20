import numpy as np
from ozy.dict_variables import common_variables,grid_variables,particle_variables,get_code_units,basic_conv
import sys
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/amr')
sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/part')

def unigrid_amr(obj, group=None, filter=None, lmax =0, n=[100,100,100], vars=['gas/density'], xmin=(0,'code_length'), xmax=(1,'code_length'), ymin=(0,'code_length'),
                 ymax=(1,'code_length'), zmin=(0,'code_length'), zmax=(1,'code_length'), angmom = np.array([0,0,1]), symlog=True):
    
    from ozy.utils import init_region,init_filter
    from ozy.group import create_new_group

    from amr2 import filtering
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
        filt = filtering.filter()
    
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