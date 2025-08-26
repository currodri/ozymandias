import os
import numpy as np
import pandas as pd
from .HaloMaker_utils import hmCatalogue,galaxyCatalogue
from hutils import py_halo_utils as phu
from tuitls import py_tree_ios as pti

class TreeMaker(object):
    def __init__(self, RunDir, BricksDir=None, HaloDir='Halos'):
        self.RunDir = RunDir
        self.BricksDir = BricksDir
        self.is_galaxies = False # The basic tree maker is for DM
        self.TreeDir = os.path.join(RunDir, 'Trees', HaloDir)
    
    def setup_TreeMaker(self,run=False,minSnapNum=0,maxSnapNum=300,snaps=None):

        # 1. Make a list of snapshots if these are not provided already
        if snaps is None:
            snaps = list(range(minSnapNum, maxSnapNum + 1))

        # 2. Create the tree directory if it does not exist
        if not os.path.exists(self.TreeDir):
            os.makedirs(self.TreeDir)

        # 3. Loop over all snapshots and their corresponding tree_bricks
        #    to make sure that we only have the timesteps when halos are found
        n = 0
        for ts in snaps:
            if self.is_galaxies:
                h = hmCatalogue(self.Obj,self.RunDir,ts.snapindex)
            else:
                h = galaxyCatalogue(self.Obj,self.RunDir,ts.snapindex)
            if os.path.isfile(h.TreeFile):
                nh, ns = phu.get_nb_halos(h.TreeFile)
                if nh > 0:
                    n += 1
        
        # 4. Write the parameter file for TreeMaker
        param_file = os.path.join(self.TreeDir,'input_TreeMaker.dat')
        f = open(param_file, 'w')
        f.write('%i  1 \n'%(n))
        for ts in snaps:
            if self.is_galaxies:
                h = hmCatalogue(self.Obj,self.RunDir,self.Obj.snapindex)
            else:
                h = galaxyCatalogue(self.Obj,self.RunDir,self.Obj.snapindex)
            if os.path.isfile(h.TreeFile):
                nh, ns = phu.get_nb_halos(h.TreeFile)
                if nh > 0:
                    f.write('\''+os.path.abspath(h.TreeFile)+'\'\n')
        f.close()

        if run:
            runCmd = OZYPATH + '/TreeMaker/program/TreeMaker ./ > runTree.log'
            here = os.getcwd()
            os.chdir(self.TreeDir)
            os.system('ulimit -s unlimited')  # avoids segfault due to large recursive routine... 
            os.system(runCmd)
            os.chdir(here)

    def read_tree(self):

        # 1. Look for file tstep_file_*. If more than one, use the one 
        #    with the most recent timestamp.
        tstep_files = [f for f in os.listdir(self.TreeDir) if f.startswith('tstep_file_')]
        if len(tstep_files) == 0:
            raise FileNotFoundError("No tstep_file found in TreeDir: %s" % self.TreeDir)
        elif len(tstep_files) > 1:
            tstep_files.sort(key=lambda x: os.path.getmtime(os.path.join(self.TreeDir, x)), reverse=True)
        tstep_file = os.path.join(self.TreeDir, tstep_files[0])
        print("Reading tree from file: %s" % tstep_file)


        # 2. Get the number of nsteps from the filename and check if it is valid.
        self.nsteps = int(tstep_file.split('_')[2].split('.')[0])
        nsteps = pti.get_nsteps(self.TreeDir, self.nsteps)
        if nsteps != self.nsteps:
            raise ValueError("Number of steps in file (%d) does not match expected (%d)" % (nsteps, self.nsteps))
        
        # 3. Start reading the tree file
        self.nids = pti.get_nids(self.TreeDir, self.nsteps, nsteps)
        self.nhalos = np.zeros(self.nsteps, dtype=np.int32)
        self.aexp = np.zeros(self.nsteps, dtype=np.float64)
        self.age_univ = np.zeros(self.nsteps, dtype=np.float64)
        pti.read_timestep_props(self.TreeDir, self.nsteps, nsteps,
                                self.nhalos, self.aexp, self.age_univ)
        nh = sum(self.nhalos)

        # 4. Get the number of the last snapshot used in the tree
        f = open("%s/input_TreeMaker.dat"%(self.TreeDir),'r')
        lines = f.readlines()
        f.close()
        last_snap = int(lines[-1].split('/')[-1].split('bricks')[-1].split("'\n")[0])  

        # 5. Read the tree structure
        IDs = np.zeros((nh,self.nids), dtype=np.int32)
        pti.read_tree_struct(self.TreeDir, self.nsteps, nsteps,
                             self.nhalos, nh, self.nids, IDs)
        

        # 6. Correct the tree timesteps since they start from the first snapshot
        #    that holds the halo.
        IDs[:,4] = IDs[:,4] + (last_snap - IDs[:,4].max())
        ID_keys = ('bush_id', 'tree_id', 'halo_id', 'halo_num', 'halo_ts',
                   'first_prog', 'next_prog', 'descendent_id', 'last_prog',
                   'host_halo_id', 'host_sub_id', 'next_sub_id')
        
        # 7. Create a DataFrame with the tree data
        pIDs = pd.DataFrame(IDs, columns=ID_keys)
        self.nprops = pti.get_nprops(self.TreeDir, self.nsteps)
        props = np.zeros((nh, self.nprops), dtype=np.float64)
        pti.read_props(self.TreeDir, self.nsteps, self.nhalos, nh, self.nprops, props)
        p_keys = ('x', 'y', 'z', 'vx', 'vy', 'vz', 'm', 'r', 'spin',
                      'rvir', 'mvir', 'tvir', 'cvel', 'dmacc', 'frag',
                      'Lx', 'Ly', 'Lz', 'ep', 'ek', 'et')
        pprops = pd.DataFrame(props, columns=p_keys)
        pIDs.set_index(pIDs.halo_id, inplace=True)
        pprops.set_index(pIDs.halo_id, inplace=True)
        self.tree = pd.concat([pIDs, pprops], axis=1)

        # 8. Add also the index with halo_num
        self.indid = pd.DataFrame(IDs[:,2:5],columns=('halo_id', 'halo_num', 'halo_ts'))
        self.indid.set_index([self.indid.halo_ts,self.indid.halo_num],inplace=True)

    def get_haloID(self,num,ramses_timestep):
        # First map the ramses snapshot number to Tree timestep
        tree_ts = self.ramsests2treets(ramses_timestep)
        return self.indid[(self.indid.halo_ts==tree_ts)&(self.indid.halo_num==int(num))].halo_id
        
    def get_main_progenitor(self, target_id):
        """Return the main progenitors of a halo with their timestep."""
        target = self.trees.loc[target_id]
        progenitors = self._get_progenitors(target_id)
        id_main = {}
        for ts in range(int(self.trees.loc[target_id].halo_ts)+1):
            prog_mass = progenitors[progenitors.halo_ts == ts].m
            if len(prog_mass):
                id_main[ts] = int(prog_mass.idxmax())                

        main_progs = pd.concat([self.trees.loc[id_main[ts]] for ts in id_main], axis=1, join='inner').T
        return main_progs

    def get_all_progenitors(self, target_id):
        """Return the main progenitors of a halo with their timestep."""

        target = self.trees.loc[target_id]
        progenitors = self._get_progenitors(target_id)

        return progenitors

    def get_target(self, hid):
        target = self.trees.loc[hid]
        return target

    def _get_progenitors(self, hid):
        target = self.trees.loc[hid]

        mask = ((self.trees.halo_id >= int(hid)) &
                (self.trees.halo_id <= int(target['last_prog'])))
        progenitors = self.trees.loc[mask].copy()
        return progenitors

    def get_full_tree(self,hid):
        ''' 
        Return full tree in which halo hid is . 
        '''
        target = self.trees.loc[hid]
        mask = (self.trees.tree_id == int(target['tree_id']))   
        fullTree = self.trees.loc[mask].copy()
        return fullTree

    def get_full_bush(self,hid):
        ''' 
        Return full bush in which halo hid is . 
        '''
        target = self.trees.loc[hid]
        mask = (self.trees.bush_id == int(target['bush_id']))   
        fullBush = self.trees.loc[mask].copy()
        return fullBush

    def make_timestep_mapping(self):
        import re

        # 1. Read the input_TreeMaker.dat file to get the mapping
        param_file = os.path.join(self.TreeDir, 'input_TreeMaker.dat')
        f = open(param_file, 'r') 
        lines = f.readlines()
        nLines = len(lines)
        self.timestep_map = np.zeros(nLines-1,dtype=int)
        for i in range(1,nLines):
            numbers = re.findall(r'\d+', lines[i])
            self.timestep_map[i-1] = int(numbers[len(numbers)-1])

        # 3. The TreeMaker ts is offset upwards, so we need to
        #    prepend the mapping array with empty elements
        last_ramses_snap = self.timestep_map[len(self.timestep_map)-1]
        offset = last_ramses_snap - len(self.timestep_map) + 1
        emptyArray = np.zeros(offset, dtype=int)
        self.timestep_map = np.insert(self.timestep_map, 0, emptyArray)

    def treets2ramsests(self,ts):
        ''' 
        Map a TreeMaker timestep to a ramses output number. This is necessary
        when there are gaps in the sequence of ramses snapshots, i.e. when
        snapshots have been deleted.
        '''
        if not hasattr(self, 'timestep_map'):
            self.make_timestep_mapping()
        return self.timestep_map[ts]

    def ramsests2treets(self,ramses_ts):
        ''' 
        Map a ramses snapshot number to a Tree timestep. This is necessary
        when there are gaps in the sequence of ramses snapshots, i.e. when
        snapshots have been deleted.
        '''
        if not hasattr(self, 'timestep_map'):
            self.make_timestep_mapping()
        for i in range(len(self.timestep_map)):
            if self.timestep_map[i] == ramses_ts:
                return i
        return -1


