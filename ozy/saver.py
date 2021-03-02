import os
import h5py
import numpy as np
from yt import YTArray, YTQuantity

blacklist = [
    'G', 'initial_mass',
    'valid', 'vel_conversion',
    'unbound_particles', '_units',
    'unit_registry_json',
    'unbound_indexes',
    'lists','dicts'
]

def _write_dataset(key, data, hd):
    hd.create_dataset(key, data=data, compression=1)

def check_and_write_dataset(obj, key, hd):
    """General function that writes a HDF5 dataset.
    """
    if not hasattr(obj, key):
        return
    if isinstance(getattr(obj, key), int):
        return
    _write_dataset(key, getattr(obj, key), hd)
    
def serialise_list(obj_list, key, hd):
    """Function that serialises a index list for objects.
    """
    if key in blacklist:
        return
    if not hasattr(obj_list[0], key):
        return
    data = _get_serialised_list(obj_list, key)
    
    _write_dataset(key, data, hd)
    
def _get_serialised_list(obj_list, key):
    tmp = []
    index = 0
    for i in obj_list:
        current_list = getattr(i, key)
        n = len(current_list)
        tmp.extend(current_list)
        setattr(i, '%s_start' % key, index)
        index += n
        setattr(i, '%s_end'   % key, index)
    return tmp

    
def serialise_attributes(obj_list, hd, hd_dicts):
    """Function that goes through a list of Group objects
    and serialises their attributes.
    """
    for k, v in obj_list[0].__dict__.items():
        if k in blacklist:
            continue
        if isinstance(v, dict):
            _write_dict(obj_list, k, v, hd_dicts)
        else:
            _write_attrib(obj_list, k, v, hd)

def _write_attrib(obj_list, k, v, hd):
    unit = False
    if isinstance(v, YTQuantity):
        data = [getattr(i, k).d for i in obj_list]
        unit = True
    elif isinstance(v, YTArray):
        # For the case of a vector
        if np.shape(v)[0] == 3:
            data = np.vstack([getattr(i, k).d for i in obj_list])
        else:
            # These are just single item arrays
            data = [getattr(i,k).d for i in obj_list]
    elif isinstance(v, np.ndarray) and np.shape(v)[0] == 3 and 'list' not in k:
        try:
            data = np.vstack([getattr(i,k) for i in obj_list])
        except:
            return
    elif isinstance(v, (int, float, bool, np.number)):
        data = [getattr(i,k) for i in obj_list]
    else:
        return
    
    _write_dataset(k, data, hd)
    if unit:
        hd[k].attrs.create('unit', str(v.units).encode('utf8'))

def _write_dict(obj_list, k, v, hd):
    for kk, vv in v.items():
        unit = False
        if isinstance(vv, (YTQuantity, YTArray)):
            data = np.array([getattr(i,k)[kk].d for i in obj_list])
            unit = True
        else:
            data = np.array([getattr(i,k)[kk] for i in obj_list])            

        _write_dataset('%s.%s' % (k,kk), data, hd)
        if unit:
            hd['%s.%s' % (k,kk)].attrs.create('unit', str(vv.units).encode('utf8'))
            
def seralise_global_attributes(obj, hd):
    """Function that goes through a OZY object and saves
    general attributes.
    """
    units = {}
    for k, v in obj.__dict__.items():
        if k in blacklist:
            continue
        if isinstance(v, (YTArray, YTQuantity)):
            hd.attrs.create(k, v.d)
            units[k] = v.units
        elif isinstance(v, str):
            hd.attrs.create(k, v.encode('utf8'))
        elif isinstance(v, (int, float, bool, np.number)):
            hd.attrs.create(k, v)
    # Now save the collected units.
    if len(units) > 0:
        uhd = hd.create_group('global_attribute_units')
        for k, v in units.items():
            uhd.attrs.create(k, str(v).encode('utf8'))

def save(obj, filename='default'):
    """This is the main function that saves a OZY catalogue to disk
    using the HDF5 format.
    """
    sim_folder = '/' + os.path.join(*obj.simulation.fullpath.split('/')[0:-1])
    snap_ID = obj.simulation.fullpath.split('/')[-1].split('_')[-1]
    
    if filename == 'default':
        filename = os.path.join(sim_folder, 'ozy_'+str(snap_ID)+'.hdf5')
        
    if os.path.isfile(filename):
        os.remove(filename)
    
    outfile = h5py.File(filename, 'w')
    outfile.attrs.create('ozy', 315)
    
    # Just to make possible to reconstruct the units used by yT.
    unit_registry = obj.yt_dataset.unit_registry.to_json()
    outfile.attrs.create('unit_registry_json', unit_registry.encode('utf8'))
    
    seralise_global_attributes(obj, outfile)
    obj.simulation._serialise(obj, outfile)
    
    if hasattr(obj, 'halos') and obj.nhalos > 0:
        hd = outfile.create_group('halo_data')
        hdd = hd.create_group('lists')
        hddd = hd.create_group('dicts')
        
        # TODO: Addapt for particle lists.
        serialise_list(obj.halos, 'galaxy_index_list', hdd)
        serialise_attributes(obj.halos, hd, hdd)
    
    if hasattr(obj, 'galaxies') and obj.ngalaxies > 0:
        hd = outfile.create_group('galaxy_data')
        hdd = hd.create_group('lists')
        hddd = hd.create_group('dicts')
        
        # TODO: Addapt for particle lists.
        index_lists = ['slist','cloud_index_list']
        for vals in index_lists:
            serialise_list(obj.galaxies, vals, hdd)
        serialise_attributes(obj.galaxies, hd, hddd)
    
    if hasattr(obj, 'clouds') and obj.nclouds > 0:
        hd = outfile.create_group('cloud_data')
        hdd = hd.create_group('lists')
        hddd = hd.create_group('dicts')
        
        # TODO: Addapt for particle lists.
        # serialise_list(obj.galaxies, 'glist', hdd)
        serialise_attributes(obj.clouds, hd, hddd)
    
    # TODO: Allow to save global_particle_lists
    # if hasattr(obj, 'global_particle_lists'):
    
    outfile.close()