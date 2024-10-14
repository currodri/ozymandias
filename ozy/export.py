import numpy as np
from ozy.dict_variables import common_variables,grid_variables,particle_variables,get_code_units,basic_conv

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

def part2chunk(obj, group=None, filter=None, outname=None, ptype='star',vars=['star/x','star/y','star/z','star/mass'],
               xmin=(0,'code_length'), xmax=(1,'code_length'), ymin=(0,'code_length'),
                 ymax=(1,'code_length'), zmin=(0,'code_length'), zmax=(1,'code_length'),rmin=(0,'kpc'),rmax=(0.5,'code_length'), angmom = np.array([0,0,1]), 
                 filetype='hdf5'):
    
    import os
    from ozy.utils import init_region
    from ozy.group import create_new_group
    
    from part2 import filtering,export_part
    
    output_path = obj.simulation.fullpath
    
    # 1. Initialise region
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
        
    # 2. Get filter if not given
    if filter == None:
        filt = filtering.filter()
        
    # 3. Set the name of the output file if not given
    outid = output_path.split('/')[-1][-5:]
    if outname == None:
        if group == None:
            outname = os.getcwd()+f'/{ptype}_'+outid+'_noneobj'
        else:
            outname = os.getcwd()+f'/{ptype}_'+outid+'_gal'+str(group.ID)
        
    # 4. Get the particle data in a chunk
    partchunk = export_part.chunk_handler()
    partchunk.nvars = len(vars)
    partchunk.filt = filt
    partchunk.ptype = ptype
    
    export_part.allocate_chunk_handler(partchunk)
    for i in range(0, len(vars)):
        partchunk.varnames.T.view('S128')[i] = vars[i].ljust(128)
        
    print('Extracting %i variables from the particles in the region...'%len(vars))
    export_part.part2chunk(output_path,selected_reg,partchunk)
    
    # 5. Save to file in the requested format
    if filetype == 'txt':
        savetotxt(partchunk,outname,vars)
    else:
        savetohdf5(partchunk,outname,vars)
        
def savetotxt(partchunk, outname, vars):
    """
    Save partchunk data to a plain text file using provided variable names.

    Parameters:
    partchunk : object
        The chunk handler object containing the data to save.
    outname : str
        The name of the output text file.
    vars : list of str
        The list of variable names to use as headers.
    """
    # Open the output text file for writing
    with open(outname + '.txt', 'w') as f:
        # Create header with variable names and units
        header_vars = ' '.join(vars)
        header_units = ' '.join([get_code_units(var.split('/')[-1]) for var in vars])
        
        # Write variable names and their units as two header lines
        f.write('# Variables: ' + header_vars + '\n')
        f.write('# Units: ' + header_units + '\n')
        
        # Loop through the particles and write each one's data
        for i in range(partchunk.npartsaved):
            # Extract data for all variables for this particle (i.e., all vars for particle i)
            line_data = ' '.join([str(partchunk.data[var_index, i]) for var_index in range(partchunk.nvars)])
            f.write(line_data + '\n')

    print(f'Data successfully saved to {outname}.txt')

def savetohdf5(partchunk, outname, vars):
    """
    Save partchunk data to an HDF5 file using provided variable names.

    Parameters:
    partchunk : object
        The chunk handler object containing the data to save.
    outname : str
        The name of the output HDF5 file.
    vars : list of str
        The list of variable names to use as dataset names.
    """
    import h5py
    # Open an HDF5 file in write mode
    with h5py.File(outname + '.hdf5', 'w') as f:
        # Loop through the variables and save each as a dataset using names from `vars`
        for var_index in range(partchunk.nvars):
            var_name = vars[var_index]
            data = partchunk.data[var_index, :partchunk.npartsaved]  # Select data for this variable
            dataset = f.create_dataset(var_name, data=data)
            
            # Add the unit as an attribute to the dataset
            unit = get_code_units(var_name.split('/')[-1])
            dataset.attrs['unit'] = unit

    print(f'Data successfully saved to {outname}.hdf5')
    

def export2skirt(obj, group=None, filter=None, lmax =0, var='dust_density', xmin=(0,'code_length'), xmax=(1,'code_length'), ymin=(0,'code_length'),
                 ymax=(1,'code_length'), zmin=(0,'code_length'), zmax=(1,'code_length'),rmin=(0,'kpc'),rmax=(0.5,'code_length'), angmom = np.array([0,0,1]), 
                 symlog=True, h=(30,'pc'), smoothmethod='constant',sedmethod='bruzual&charlot',recompute=False):
    
    import os
    import ozy
    from ozy.utils import init_region,init_filter
    from ozy.group import create_new_group

    from part2 import filtering,export_part
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
        filt = filtering.filter()
    
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
    from ozy.utils import init_region,init_filter
    from ozy.group import create_new_group

    from part2 import filtering,export_part

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
        filt = filtering.filter()
    
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
