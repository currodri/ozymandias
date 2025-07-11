import functools
import os.path
from pprint import pprint

import h5py
import numpy as np
from unyt import UnitRegistry,unyt_array,unyt_quantity
from collections import defaultdict
from .sim_attributes import SimulationAttributes
from .utils import info_printer
from .saver import _write_attrib, _write_dict
from .lazy_classes import LazyList,LazyDataset,LazyDict
from unyt.dimensions import length,mass,time,temperature
from unyt import mp,erg,K,g,kb


class Profile:
    def __init__(self, obj, index, group_type, key, hd, nosubs):
        self.obj = obj
        self._index = index
        self.group_type = group_type
        self.key = key
        self.nbins = 0
        self.xvar = None
        self.region = None
        self.filter = {}
        self.lmax = 0
        self.yvars = {}
        self.weightvars = {}
        self.xdata = None
        self.ydata = {}
        self.nosubs = nosubs

        self.blacklist = [
            'obj','_index','group_type','region',
            'filter','yvars','weightvars','xdata','ydata'
        ]
        self._unpack(hd)
    
    def _unpack(self, hd):
        if self.nosubs:
            path = str(self.group_type+'_data/profiles_nosubs/'+str(self._index)+'/'+self.key)
        else:
            path = str(self.group_type+'_data/profiles/'+str(self._index)+'/'+self.key)
        profile_gp = hd[path]
        for k,v in profile_gp.attrs.items():
            if k in self.blacklist:
                continue
            else:
                setattr(self, k, v)
        self.region = {}
        self.region['type'] = profile_gp.attrs['type']
        self.region['centre'] = profile_gp.attrs['centre']
        self.region['axis'] = profile_gp.attrs['axis']
        self.filter = {}
        self.filter['name'] = profile_gp.attrs['name']
        self.filter['conditions'] = profile_gp['conditions'][:]

        for k in profile_gp.keys():
            if k != 'xdata' and k != 'conditions':
                self.yvars[k] = []
                self.ydata[k] = []
                self.weightvars[k] = []
                for j in profile_gp[k].keys():
                    self.yvars[k].append(j)
                    data = profile_gp[k+'/'+j]
                    self.ydata[k].append(unyt_array(data[:], str(data.attrs['units']), registry=self.obj.unit_registry))
                    self.weightvars[k] = list(data.attrs['weightvars'])
            else:
                self.xdata = unyt_array(profile_gp['xdata'][:], profile_gp['xdata'].attrs['units'], registry=self.obj.unit_registry)

class PhaseDiagram:
    def __init__(self,obj,index,group_type,key,hd,nosubs):
        self.obj = obj
        self._index = index
        self.group_type = group_type
        self.key = key
        self.nbins = [0,0]
        self.xvar = None
        self.yvar = None
        self.region = None
        self.filter = {}
        self.lmax = 0
        self.zvars = {}
        self.weightvars = {}
        self.xdata = None
        self.ydata = None
        self.zdata = {}
        self.nosubs = nosubs

        self.blacklist = [
            'obj','_index','group_type','region',
            'filter','zvars','weightvars','xdata','ydata'
        ]
        self._unpack(hd)

    def _unpack(self, hd):
        if self.nosubs:
            path = str(self.group_type+'_data/histograms_nosubs/'+str(self._index)+'/'+self.key)
        else:
            path = str(self.group_type+'_data/histograms/'+str(self._index)+'/'+self.key)            
        histograms_gp = hd[path]
        for k,v in histograms_gp.attrs.items():
            if k in self.blacklist:
                continue
            else:
                setattr(self, k, v)
        self.region = {}
        self.region['type'] = histograms_gp.attrs['type']
        self.region['centre'] = histograms_gp.attrs['centre']
        self.region['axis'] = histograms_gp.attrs['axis']
        self.filter = {}
        self.filter['name'] = histograms_gp.attrs['name']
        self.filter['conditions'] = histograms_gp['conditions'][:]
        for k in histograms_gp.keys():
            if k != 'xdata' and k != 'ydata' and k != 'conditions':
                self.zvars[k] = []
                self.zdata[k] = []
                for j in histograms_gp[k].keys():
                    if j == 'weightvars':
                        self.weightvars[k] = histograms_gp[k+'/'+j][:]
                        for w in range(0, len(self.weightvars[k])):
                            self.weightvars[k][w] = self.weightvars[k][w].decode('utf8')
                        continue
                    self.zvars[k].append(j)
                    data = histograms_gp[k+'/'+j]
                    self.zdata[k].append(unyt_array(data[:], str(data.attrs['units']), registry=self.obj.unit_registry))
            else:
                self.xdata = unyt_array(histograms_gp['xdata'][:], histograms_gp['xdata'].attrs['units'], registry=self.obj.unit_registry)
                self.ydata = unyt_array(histograms_gp['ydata'][:], histograms_gp['ydata'].attrs['units'], registry=self.obj.unit_registry)

class GalacticFlow:
    def __init__(self,obj,index,group_type,key,hd,nosubs):
        self.obj = obj
        self._index = index
        self.group_type = group_type
        self.key = key
        self.region = None
        self.type = 'none'
        self.filter = {}
        self.data = {}
        self.nosubs = nosubs

        self.blacklist = [
            'obj','_index','group_type','region',
            'filter'
        ]
        self._unpack(hd)

    def _unpack(self,hd):
        if self.nosubs:
            path = str(self.group_type+'_data/flows_nosubs/'+str(self._index)+'/'+self.key)
        else:
            path = str(self.group_type+'_data/flows/'+str(self._index)+'/'+self.key)
        flows_gp = hd[path]
        for k,v in flows_gp.attrs.items():
            if k not in self.blacklist:
                setattr(self, k, v)
        self.region = {}
        self.region['type'] = flows_gp.attrs['type']
        self.region['centre'] = flows_gp.attrs['centre']
        self.region['axis'] = flows_gp.attrs['axis']
        self.region['r'] = unyt_quantity(flows_gp.attrs['r'], 'code_length', registry=self.obj.unit_registry)
        self.region['rmin'] = unyt_quantity(flows_gp.attrs['rmin'], 'code_length', registry=self.obj.unit_registry)
        self.region['rmax'] = unyt_quantity(flows_gp.attrs['rmax'], 'code_length', registry=self.obj.unit_registry)
        self.region['dr'] = unyt_quantity(flows_gp.attrs['dr'], 'code_length', registry=self.obj.unit_registry)
        self.filter = {}
        self.type = flows_gp.attrs['name']
        self.filter['name'] = flows_gp.attrs['name']
        self.filter['conditions'] = flows_gp['conditions'][:]

        for k in flows_gp.keys():
            if k != 'conditions':
                data = flows_gp[k]
                for j in data.keys():
                    unit = data[j].attrs['units']
                    if data[j].shape == 1:
                        if unit =='code_energy':
                            unit = 'code_mass * code_velocity**2'
                        self.data[j] = unyt_quantity(data[j][()], unit, registry=self.obj.unit_registry)
                    else:
                        if unit =='code_energy':
                            unit = 'code_mass * code_velocity**2'
                        self.data[j] = unyt_array(data[j][()], unit, registry=self.obj.unit_registry)

class OZY:
    def __init__(self, filename, read_mode='r'):
        self._ds = None
        self.data_file = os.path.abspath(filename)
        self.snapname = self.data_file.split('/')[-1]
        self.snapID = self.snapname.split('.')[0].split('_')[-1]
        self._set_variable_ordering()
        self._galaxy_slist = LazyDataset(self, 'galaxy_data/lists/slist')

        hd = h5py.File(filename, read_mode)
        # This should be the ozy_version with which the dataset was created.
        self.ozy = hd.attrs['ozy']
        self.unit_registry = UnitRegistry.from_json(
            hd.attrs['unit_registry_json'])

        # Load the simulation attributes.
        self.simulation = SimulationAttributes()
        self.simulation._unpack(self, hd)
        # self.simulation._update(self, hd)
        
        
        # Halo data is loaded ALWAYS.
        # TODO: Change this so that it's just a flag.
        self._galaxy_index_list = None
        if 'halo_data/lists/galaxy_index_list' in hd:
            self._galaxy_index_list = LazyDataset(
                self, 'halo_data/lists/galaxy_index_list')
        
        self._halo_data = {}
        for k, v in hd['halo_data'].items():
            if type(v) is h5py.Dataset:
                self._halo_data[k] = LazyDataset(self, 'halo_data/' + k)

        self._halo_dicts = defaultdict(dict)
        for k in hd['halo_data/dicts']:
            dictname, arrname = k.split('.')
            self._halo_dicts[dictname][arrname] = LazyDataset(
                self, 'halo_data/dicts/' + k)
        
        self.nhalos = hd.attrs['nhalos']
        self.halos  = LazyList(self.nhalos, lambda i: Halo(self, i))
        
        # Not all snaps will have galaxies, so we need to load firstly 
        # default values for everything.
        self.have_flows = False
        self.have_flows_nosubs = False
        self.have_profiles = False
        self.have_profiles_nosubs = False
        self.have_histograms = False
        self.have_histograms_nosubs = False
        self._galaxy_data  = {}
        self._galaxy_dicts = defaultdict(dict)
        self.ngalaxies     = 0
        self.galaxies      = LazyList(self.ngalaxies, lambda i: Galaxy(self, i))
        if 'galaxy_data' in hd:
            self._cloud_index_list = None
            if 'galaxy_data/lists/cloud_index_list' in hd:
                self._cloud_index_list = LazyDataset(
                    self, 'galaxy_data/lists/cloud_index_list')

            if 'tree_data/progen_galaxy_star' in hd:
                self._galaxy_data['progen_galaxy_star'] = self._progen_galaxy_star = LazyDataset(
                    self, 'tree_data/progen_galaxy_star')
                
            if 'tree_data/descend_galaxy_star' in hd:
                self._galaxy_data['descend_galaxy_star'] = self._descend_galaxy_star = LazyDataset(
                    self, 'tree_data/descend_galaxy_star')

            for k, v in hd['galaxy_data'].items():
                if type(v) is h5py.Dataset:
                    self._galaxy_data[k] = LazyDataset(
                        self, 'galaxy_data/' + k)

            for k in hd['galaxy_data/dicts']:
                dictname, arrname = k.split('.')
                self._galaxy_dicts[dictname][arrname] = LazyDataset(
                    self, 'galaxy_data/dicts/' + k)

            self.ngalaxies = hd.attrs['ngalaxies']
            self.galaxies = LazyList(self.ngalaxies,
                                        lambda i: Galaxy(self, i))

            if 'galaxy_data/profiles' in hd:
                self.have_profiles = True
                prof_indices = []
                prof_keys = []
                for k in hd['galaxy_data/profiles'].keys():
                    for j in hd['galaxy_data/profiles/'+k].keys():
                        prof_indices.append(int(k))
                        prof_keys.append(j)
                self._galaxy_profile_index_list = LazyList(
                    len(prof_indices), lambda i: int(prof_indices[i])
                )
                self._galaxy_profiles = [Profile(self,int(prof_indices[i]),'galaxy', prof_keys[i],hd,False) for i in range(0, len(prof_indices))]
            if 'galaxy_data/profiles_nosubs' in hd:
                self.have_profiles_nosubs = True
                prof_nosubs_indices = []
                prof_nosubs_keys = []
                for k in hd['galaxy_data/profiles_nosubs'].keys():
                    for j in hd['galaxy_data/profiles_nosubs/'+k].keys():
                        prof_nosubs_indices.append(int(k))
                        prof_nosubs_keys.append(j)
                self._galaxy_profile_nosubs_index_list = LazyList(
                    len(prof_nosubs_indices), lambda i: int(prof_nosubs_indices[i])
                )
                self._galaxy_profiles_nosubs = [Profile(self,int(prof_nosubs_indices[i]),'galaxy', prof_nosubs_keys[i],hd,True) for i in range(0, len(prof_nosubs_indices))]
            if 'galaxy_data/histograms' in hd:
                self.have_histograms = True
                pd_indices = []
                pd_keys = []
                for k in hd['galaxy_data/histograms'].keys():
                    for j in hd['galaxy_data/histograms/'+k].keys():
                        pd_indices.append(int(k))
                        pd_keys.append(j)
                self._galaxy_phasediag_index_list = LazyList(
                    len(pd_indices), lambda i: int(pd_indices[i])
                )
                self._galaxy_phasediag = [PhaseDiagram(self,int(pd_indices[i]),'galaxy', pd_keys[i],hd,False) for i in range(0, len(pd_indices))]
            if 'galaxy_data/histograms_nosubs' in hd:
                self.have_histograms_nosubs = True
                pd_nosubs_indices = []
                pd_nosubs_keys = []
                for k in hd['galaxy_data/histograms_nosubs'].keys():
                    for j in hd['galaxy_data/histograms_nosubs/'+k].keys():
                        pd_nosubs_indices.append(int(k))
                        pd_nosubs_keys.append(j)
                self._galaxy_phasediag_nosubs_index_list = LazyList(
                    len(pd_nosubs_indices), lambda i: int(pd_nosubs_indices[i])
                )
                self._galaxy_phasediag_nosubs = [PhaseDiagram(self,int(pd_nosubs_indices[i]),'galaxy', pd_nosubs_keys[i],hd,True) for i in range(0, len(pd_nosubs_indices))]
            if 'galaxy_data/flows' in hd:
                self.have_flows = True
                gf_indices = []
                gf_keys = []
                for k in hd['galaxy_data/flows'].keys():
                    for j in hd['galaxy_data/flows/'+k].keys():
                        gf_indices.append(int(k))
                        gf_keys.append(j)
                self._galaxy_flows_index_list = LazyList(
                    len(gf_indices), lambda i: int(gf_indices[i])
                )
                self._galaxy_flows = [GalacticFlow(self,int(gf_indices[i]),'galaxy',gf_keys[i],hd,False) for i in range(0, len(gf_indices))]
            if 'galaxy_data/flows_nosubs' in hd:
                self.have_flows_nosubs = True
                gf_nosubs_indices = []
                gf_nosubs_keys = []
                for k in hd['galaxy_data/flows_nosubs'].keys():
                    for j in hd['galaxy_data/flows_nosubs/'+k].keys():
                        gf_nosubs_indices.append(int(k))
                        gf_nosubs_keys.append(j)
                self._galaxy_flows_nosubs_index_list = LazyList(
                    len(gf_nosubs_indices), lambda i: int(gf_nosubs_indices[i])
                )
                self._galaxy_flows_nosubs = [GalacticFlow(self,int(gf_nosubs_indices[i]),'galaxy',gf_nosubs_keys[i],hd,True) for i in range(0, len(gf_nosubs_indices))]
    
        hd.close()
        del hd
    @property
    def central_galaxies(self):
        return [h.central_galaxy for h in self.halos]

    @property
    def satellite_galaxies(self):
        galaxies = []
        for h in self.halos:
            galaxies.extend(h.satellite_galaxies)
    
    def most_massive_galaxy(self,return_index=False):
        mstellar = [i.mass['stellar'] for i in self.galaxies]
        ibig = np.argmax(mstellar)
        if return_index:
            return ibig
        else:
            return self.galaxies[ibig]
        
    def most_massive_system(self,return_index=False):
        mvir = [i.virial_quantities['mass'] for i in self.galaxies]
        ibig = np.argmax(mvir)
        if return_index:
            return ibig
        else:
            return self.galaxies[ibig]
    def array(self, value, units):
        return unyt_array(value, units, registry=self.unit_registry)

    def quantity(self, value, units):
        return unyt_quantity(value, units, registry=self.unit_registry)

    def galinfo(self, top=10):
        info_printer(self, 'galaxy', top)

    def haloinfo(self, top=10):
        info_printer(self, 'halo', top)

    def cloudinfo(self, top=10):
        info_printer(self, 'cloud', top)

    def _set_variable_ordering(self):
        from amr2 import dictionary_commons
        from variables_settings import hydro_variables_ordering, \
                                        part_variables_ordering, \
                                        part_variables_type
        self.vardict = dictionary_commons.dictf90()
        self.part_vardict = dictionary_commons.dictf90()
        self.part_vartypes = dictionary_commons.dictf90()
        if len(hydro_variables_ordering) > 0:
            print("Setting hydro variable ordering from local ozy_settings.py.")
            self.vardict.init(len(hydro_variables_ordering))
            for varname,varindex in hydro_variables_ordering.items():
                self.vardict.add(varname,varindex)
            self.use_vardict = True
        else:
            self.use_vardict = False
        if len(part_variables_ordering) > 0:
            print("Setting particle variable ordering from local ozy_settings.py.")
            self.part_vardict.init(len(part_variables_ordering))
            for varname,varindex in part_variables_ordering.items():
                self.part_vardict.add(varname,varindex)
            self.use_part_vardict = True
            self.part_vartypes.init(len(part_variables_type))
            for var,vtype in part_variables_type.items():
                self.part_vartypes.add(var,vtype)
        else:
            self.use_part_vardict = False

# TODO: All of these classes need to be adapted to be able to load particle lists.
        
class Group:
    
    def info(self):
        pdict = {}
        for k in getattr(self.obj, '_{}_data'.format(self.type)):
            pdict[k] = getattr(self, k)
        for k in getattr(self.obj, '_{}_dicts'.format(self.type)):
            pdict[k] = dict(getattr(self, k))
        pprint(pdict)
        
class Halo(Group):
    def __init__(self, obj, index):
        self.obj_type = 'halo'
        self.obj = obj
        self._index = index
        self._galaxies = None
        self._satellite_galaxies = None
        self._central_galaxy = None

    def __dir__(self):
        return dir(type(self)) + list(self.__dict__) + list(
            self.obj._halo_data) + list(self.obj._halo_dicts)
        
    @property
    def galaxy_index_list(self):
        return self.obj._galaxy_index_list[self.galaxy_index_list_start:self.
                                           galaxy_index_list_end]

    def _init_galaxies(self):
        self._galaxies = []
        self._satellite_galaxies = []
        for galaxy_index in self.galaxy_index_list:
            galaxy = self.obj.galaxies[galaxy_index]
            self._galaxies.append(galaxy)
            if galaxy.central:
                self._central_galaxy = galaxy
            else:
                self._satellite_galaxies.append(galaxy)

    @property
    def galaxies(self):
        if self._galaxies is None:
            self._init_galaxies()
        return self._galaxies

    @property
    def central_galaxy(self):
        if self._central_galaxy is None:
            self._init_galaxies()
        return self._central_galaxy

    @property
    def satellite_galaxies(self):
        if self._satellite_galaxies is None:
            self._init_galaxies()
        return self._satellite_galaxies

    @property
    def substructure_list(self):
        subs = []
        if self.nextsub == 0:
            print('This halo does not seem to have a substructure assigned!')
            return subs
        haloIDs = self.obj._halo_data['ID'][:]
        nexti = self.nextsub
        while nexti != -1:
            try:
                nextindex = np.where(nexti == haloIDs)[0][0]
                subs.append(self.obj.halos[nextindex])
                nexti = self.obj.halos[nextindex].nextsub
            except:
                nexti = -1
        return subs


    @functools.lru_cache(maxsize=None)
    def __getattr__(self, attr):
        if attr in self.obj._halo_data:
            return self.obj._halo_data[attr][self._index]
        if attr in self.obj._halo_dicts:
            return LazyDict(
                self.obj._halo_dicts[attr].keys(),
                lambda d: self.obj._halo_dicts[attr][d][self._index])
        raise AttributeError("'{}' object has no attribute '{}'".format(
            self.__class__.__name__, attr))

class Galaxy(Group):
    def __init__(self, obj, index):
        self.type = 'galaxy'
        self.obj = obj
        self._index = index
        self.halo = obj.halos[self.parent_halo_index]
        self.profiles = None
        self.histograms = None
        self.flows = None

    def __dir__(self):
        return dir(type(self)) + list(self.__dict__) + list(
            self.obj._galaxy_data) + list(self.obj._galaxy_dicts)

    def _init_profiles(self):
        self.profiles = []
        if self.obj.have_profiles:
            for p,profile_index in enumerate(self.obj._galaxy_profile_index_list):
                if profile_index == self._index:
                    profile = self.obj._galaxy_profiles[p]
                    self.profiles.append(profile)
                
        self.profiles_nosubs = []
        if self.obj.have_profiles_nosubs:
            for p,profile_index in enumerate(self.obj._galaxy_profile_nosubs_index_list):
                if profile_index == self._index:
                    profile = self.obj._galaxy_profiles_nosubs[p]
                    self.profiles_nosubs.append(profile)
    
    def _init_histograms(self):
        self.histograms = []
        if self.obj.have_histograms:
            for p,phasediag_index in enumerate(self.obj._galaxy_phasediag_index_list):
                if phasediag_index == self._index:
                    phasediag = self.obj._galaxy_phasediag[p]
                    self.histograms.append(phasediag)
        
        self.histograms_nosubs = []
        if self.obj.have_histograms_nosubs:
            for p,phasediag_index in enumerate(self.obj._galaxy_phasediag_nosubs_index_list):
                if phasediag_index == self._index:
                    phasediag = self.obj._galaxy_phasediag_nosubs[p]
                    self.histograms_nosubs.append(phasediag)

    def _init_flows(self):
        self.flows = []
        if self.obj.have_flows:
            for g,flow_index in enumerate(self.obj._galaxy_flows_index_list):
                if flow_index == self._index:
                    flow = self.obj._galaxy_flows[g]
                    self.flows.append(flow)

        self.flows_nosubs = []
        if self.obj.have_flows_nosubs:
            for g,flow_index in enumerate(self.obj._galaxy_flows_nosubs_index_list):
                if flow_index == self._index:
                    flow = self.obj._galaxy_flows_nosubs[g]
                    self.flows_nosubs.append(flow)

    @property
    def slist(self):
        return self.obj._galaxy_slist[self.slist_start:self.slist_end]

    @property
    def satellites(self):
        if self.central:
            return self.halo.satellite_galaxies
        return []

    @functools.lru_cache(maxsize=None)
    def __getattr__(self, attr):
        if attr in self.obj._galaxy_data:
            return self.obj._galaxy_data[attr][self._index]
        if attr in self.obj._galaxy_dicts:
            return LazyDict(
                self.obj._galaxy_dicts[attr].keys(),
                lambda d: self.obj._galaxy_dicts[attr][d][self._index])
        raise AttributeError("'{}' object has no attribute '{}'".format(
            self.__class__.__name__, attr))
        
    def clear_cache(self):
        """Clear the cache for the __getattr__ method."""
        self.__getattr__.cache_clear()
        
    def delete_attribute(self,attr):
        if attr in self.obj._galaxy_data:
            # A. The attribute is a single data point (no dict)
            del self.obj._galaxy_data[attr]
            hd = h5py.File(self.obj.data_file, 'r+')
            del hd[f'galaxy_data/{attr}']
            hd.close()
            self.clear_cache()
        elif attr in self.obj._galaxy_dicts:
            # B. The attribute is a dictionary with keys and values
            hd = h5py.File(self.obj.data_file, 'r+')
            for kk in self.obj._galaxy_dicts[attr].keys():
                del hd[f'galaxy_data/dicts/{attr}.{kk}']
            hd.close()
            del self.obj._galaxy_dicts[attr]
            self.clear_cache()
        elif attr.split('.')[0] in self.obj._galaxy_dicts \
            and attr.split('.')[1] in self.obj._galaxy_dicts[attr.split('.')[0]]:
            # C. The attribute is a single key of a dictionary already present
            dictname = attr.split('.')[0]
            keyname = attr.split('.')[1]
            hd = h5py.File(self.obj.data_file, 'r+')
            del hd[f'galaxy_data/dicts/{dictname}.{keyname}']
            hd.close()
            del self.obj._galaxy_dicts[dictname][keyname]
            self.clear_cache()
        
    def update_attribute(self,attr,value):
        if attr in self.obj._galaxy_data:
            # A. The attribute is a single data point (no dict)
            # 1. We update the value in the currently loaded OZY object
            data = self.obj._galaxy_data[attr][:]
            data[self._index] = value
            self.clear_cache()
            
            # 2. Update the original HDF5 file for future reference
            with h5py.File(self.obj.data_file, 'a') as hd:
                _write_attrib(self.obj.galaxies, attr, value, hd['galaxy_data'])
        elif attr in self.obj._galaxy_dicts and isinstance(value, dict):
            # B. The attribute is a dictionary with keys and values
            # 1. Update the full dictionary in the currently loaded OZY object
            for kk, vv in value.items():
                data = self.obj._galaxy_dicts[attr][kk][:]
                data[self._index] = vv
            self.clear_cache()
            
            # 2. Update the original HDF5 file for future reference
            with h5py.File(self.obj.data_file, 'a') as hd:
                _write_dict(self.obj.galaxies, attr, value, hd['galaxy_data/dicts'])
        elif attr.split('.')[0] in self.obj._galaxy_dicts \
            and attr.split('.')[1] in self.obj._galaxy_dicts[attr.split('.')[0]]:
            # C. The attribute is a single key of a dictionary already present
            # 1. Update the single key of the dictionary in the currently loaded OZY object
            dictname = attr.split('.')[0]
            keyname = attr.split('.')[1]
            data = self.obj._galaxy_dicts[dictname][keyname][:]
            data[self._index] = value
            self.clear_cache()
            # 2. Update the original HDF5 file for future reference
            with h5py.File(self.obj.data_file, 'a') as hd:
                _write_dict(self.obj.galaxies, dictname, self.obj._galaxy_dicts[dictname], hd['galaxy_data/dicts'])
        else:
            # D. This is the case of an inexistent attribute/dict
            #    Choose what type of data is it
            if isinstance(value, dict):
                # This is the case of adding a new dictionary
                # 1. Add the new data to the global OZY object
                for kk, vv in value.items():
                    k = attr + '.' + kk
                    self.obj._galaxy_dicts[attr][kk] = LazyDataset(
                        self.obj, 'galaxy_data/dicts/' + k
                    )
                    if isinstance(vv, unyt_quantity):
                        empty_array = np.full(self.obj.ngalaxies, 0.0)
                        self.obj._galaxy_dicts[attr][kk]._data = self.obj.array(empty_array, vv.units)
                    elif isinstance(vv, unyt_array):
                        empty_array = np.full((self.obj.ngalaxies,)+np.shape(vv), 0.0)
                        self.obj._galaxy_dicts[attr][kk]._data = self.obj.array(empty_array, vv.units)
                    elif isinstance(vv, (np.ndarray,list)):
                        empty_array = np.full((self.obj.ngalaxies,)+np.shape(vv), 0.0)
                        self.obj._galaxy_dicts[attr][kk]._data = empty_array
                    else:
                        empty_array = np.zeros(self.obj.ngalaxies, dtype=type(value))
                        self.obj._galaxy_dicts[attr][kk]._data = empty_array
                    self.obj._galaxy_dicts[attr][kk][self._index] = vv
                    
                # 2. Clean cache to reload Galaxy details
                self.clear_cache()
                
                # 3. Update the original HDF5 file
                with h5py.File(self.obj.data_file, 'r+') as hd:
                    _write_dict(self.obj.galaxies, attr, value, hd['galaxy_data/dicts'])
            elif len(attr.split('.')) > 1:
                # This is the case that we want to add a key to an
                # already existent dictionary
                dictname = attr.split('.')[0]
                keyname = attr.split('.')[1]
                # 1. Add the new data to the global OZY object
                self.obj._galaxy_dicts[dictname][keyname] = LazyDataset(
                        self.obj, 'galaxy_data/dicts/' + attr
                    )
                if isinstance(value, unyt_quantity):
                    empty_array = np.full(self.obj.ngalaxies, 0.0)
                    self.obj._galaxy_dicts[dictname][keyname]._data = self.obj.array(empty_array, value.units)
                    self.obj._galaxy_dicts[dictname][keyname].units = value.units
                elif isinstance(value, unyt_array):
                    empty_array = np.full((self.obj.ngalaxies,)+np.shape(value), 0.0)
                    self.obj._galaxy_dicts[dictname][keyname]._data = self.obj.array(empty_array, value.units)
                    self.obj._galaxy_dicts[dictname][keyname].units = value.units
                elif isinstance(value, (np.ndarray,list)):
                    empty_array = np.full((self.obj.ngalaxies,)+np.shape(value), 0.0)
                    self.obj._galaxy_dicts[dictname][keyname]._data = empty_array
                else:
                    empty_array = np.zeros(self.obj.ngalaxies, dtype=type(value))
                    self.obj._galaxy_dicts[dictname][keyname]._data = empty_array
                self.obj._galaxy_dicts[dictname][keyname][self._index] = value

                # 2. Clean cache to reload Galaxy details
                self.clear_cache()
                
                # 3. Update the original HDF5 file
                with h5py.File(self.obj.data_file, 'r+') as hd:
                    _write_dict(self.obj.galaxies, dictname, self.obj._galaxy_dicts[dictname], hd['galaxy_data/dicts'])
                
            else:
                # This is the final case in which we only want to add
                # a single attribute
                # 1. Add the new data to the global OZY object
                self.obj._galaxy_data[attr] = LazyDataset(
                        self.obj, 'galaxy_data/' + attr
                    )
                if isinstance(value, unyt_quantity):
                    empty_array = np.full(self.obj.ngalaxies, 0.0)
                    self.obj._galaxy_data[attr]._data = self.obj.array(empty_array, value.units)
                elif isinstance(value, unyt_array):
                    empty_array = np.full((self.obj.ngalaxies,)+np.shape(value), 0.0)
                    self.obj._galaxy_data[attr]._data = self.obj.array(empty_array, value.units)
                elif isinstance(value, (np.ndarray,list)):
                    empty_array = np.full((self.obj.ngalaxies,)+np.shape(value), 0.0)
                    self.obj._galaxy_data[attr]._data = empty_array
                else:
                    empty_array = np.zeros(self.obj.ngalaxies, dtype=type(value))
                    self.obj._galaxy_data[attr]._data = empty_array
                self.obj._galaxy_data[attr][self._index] = value

                # 2. Clean cache to reload Galaxy details
                self.clear_cache()
                
                with h5py.File(self.obj.data_file, 'r+') as hd:
                    _write_attrib(self.obj.galaxies, attr, value, hd['galaxy_data'])
            

class Cloud(Group):
    def __init__(self, obj, index):
        self.type = 'cloud'
        self.obj = obj
        self._index = index
        self.galaxy = obj.galaxies[self.parent_galaxy_index]
        self.halo = self.galaxy.halo

    def __dir__(self):
        return dir(type(self)) + list(self.__dict__) + list(
            self.obj._cloud_data) + list(self.obj._cloud_dicts)

    @functools.lru_cache(maxsize=None)
    def __getattr__(self, attr):
        if attr in self.obj._cloud_data:
            return self.obj._cloud_data[attr][self._index]
        if attr in self.obj._cloud_dicts:
            return LazyDict(
                self.obj._cloud_dicts[attr].keys(),
                lambda d: self.obj._cloud_dicts[attr][d][self._index])
        raise AttributeError("'{}' object has no attribute '{}'".format(
            self.__class__.__name__, attr))
        

# FINALLY, the function that we want!
def load(filename):
    return OZY(filename)
