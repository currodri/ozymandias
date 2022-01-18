import os
import h5py

import ozy
from ozy.progen import progen_build

from ozy.read_VELOCIraptor import make_forest, make_walkable_tree, requestedfields, additionalrequestedfields

# Load VELOCIraptor python routines (written by the developers of VELOCIraptor)
# Load the cythonized code if compiled
try:
    import velociraptor_python_tools_cython as vpt
    print('Using cython VR+TF toolkit')
except:
    import velociraptor_python_tools as vpt
    print('Using python VR+TF toolkit')

class Snapshot(object):
    """Class for tracking paths and data for simulation snapshots.
    """
    def __init__(self, snapdir, snapname, snapindex):
        self.snapdir   = snapdir
        self.snapname  = snapname
        self.snapindex = snapindex
        self.snap      = '%s/%s%05d' % (snapdir, snapname, snapindex)
        self.snapinfo  = '%s/%s%05d/info_%05d.txt' % (snapdir, snapname, snapindex, snapindex) # Assuming RAMSES convention of integer of length 5 for the snapindex
    
    def set_output_information(self, prefix='ozy_', extension='hdf5',**kwargs):
        
        self.outdir = '%s/Groups' % ('/' + os.path.join(*self.snap.split('/')[0:-1]))
        
        self.outfile = '%s/%s%s%05d.%s' % (self.outdir, prefix, self.snapname.replace(self.snapname,''), self.snapindex, extension)
    
    def _make_output_dir(self):
        if not os.path.isdir(self.outdir):
            try:
                os.makedirs(self.outdir)
            except:
                pass
    def build_HaloMaker(self, skipran, **kwargs):
        if not os.path.exists(self.snap):
            return
        self.set_output_information(**kwargs)
        
        if os.path.isfile(self.outfile) and skipran:
            return
        self._make_output_dir()
        
        obj = ozy.OZY(self.snap)
        obj.build_HaloMaker(**kwargs)
        obj.save(self.outfile)
        
        obj = None
    def build_STF(self, skipran, dm_forestdata, stars_forestdata, ind, **kwargs):
        if not os.path.exists(self.snap):
            return
        self.set_output_information(**kwargs)
        
        if os.path.isfile(self.outfile) and skipran:
            return
        self._make_output_dir()

        obj = ozy.OZY(self.snap)
        obj.build_STF(dm_forestdata, stars_forestdata, ind, **kwargs)
        obj.save(self.outfile)
        
def print_art():
    from art import text2art

    from ozy.__version__ import VERSION
    copywrite = '   (C) 2021 F. Rodriguez Montero'
    version   = '   Version %s' % VERSION

    art =  text2art("Ozymandias","ogre")
    print('\n%s\n%s\n%s\n' % (art, copywrite, version))

def drive(snapdirs, snapname, snapindexes, progen=False, skipran=False,
          halofinder='VELOCIraptor', extension='hdf5', prefix='ozy_', **kwargs):
    """Driver function for running `ÒZYMANDIAS``on multiple snapshots.
    
    """
    # Just in the case we are given a single simulation directory or
    # a single snap…
    if isinstance(snapdirs, str):
        snapdirs = [snapdirs]
    if isinstance(snapindexes, int):
        snapindexes = [int]
    
    using_mpi = False
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        nprocs = comm.Get_size()
        rank = comm.Get_rank()
        using_mpi = True
    except:
        nprocs = 1
        rank = 0
    
    if rank == 0:
        print_art()
    
    if halofinder == 'HaloMaker':
        snaps = []
        for snapdir in snapdirs:
            for snapindex in snapindexes:
                snaps.append(Snapshot(snapdir, snapname, snapindex))
        rank_snaps = snaps[rank::nprocs]
        for snap in rank_snaps:
            snap.build_HaloMaker(skipran, **kwargs)
        
        if progen:
            ozy.progen.run_progen(snapdirs, snapname, snapindexes, prefix=prefix, extension=extension, **kwargs)
    elif halofinder == 'VELOCIraptor':
        # Using a single process, check for a forest file that includes all snapshots of interest
        for snapdir in snapdirs:
            snaps = []
            for snapindex in snapindexes:
                snaps.append(Snapshot(snapdir, snapname, snapindex))
            rank_snaps = snaps[rank::nprocs]
            if rank == 0:
                stf_folder = snapdir + '/stfcat/'
                # If it's missing, create walkable tree and forest files
                if not os.path.exists(stf_folder + 'dm_forest.hdf5.0'):
                    print('DM forest file not existent. Now computing...')
                    make_walkable_tree(stf_folder,'dm')
                    make_forest(stf_folder,True,'dm')
                if not os.path.exists(stf_folder + 'stars_forest.hdf5.0'):
                    print('Stars forest file not existent. Now computing...')
                    make_walkable_tree(stf_folder,'stars')
                    make_forest(stf_folder,True,'stars')
                dm_forestdata = h5py.File(stf_folder + 'dm_forest.hdf5.0')
                stars_forestdata = h5py.File(stf_folder + 'stars_forest.hdf5.0')
                nsnaps_dm = dm_forestdata['Header'].attrs['NSnaps']
                nsnaps_stars = stars_forestdata['Header'].attrs['NSnaps']
                dm_forestdata.close()
                stars_forestdata.close()
                if nsnaps_dm != len(snapindexes):
                    print('DM forest file missing snapshots. Recomputing...')
                    make_walkable_tree(stf_folder,'dm')
                    make_forest(stf_folder,True,'dm')
                if nsnaps_stars != len(snapindexes):
                    print('Stars forest file missing snapshots. Recomputing...')
                    make_walkable_tree(stf_folder,'stars')
                    make_forest(stf_folder,True,'stars')
                
                # Assign halo and tree data to each ozy file
                fields = requestedfields + additionalrequestedfields
                # halodata, numhalos, atime, simdata, unitdata, snapnames
                dm_forestdata = vpt.ReadForest(stf_folder + 'dm_forest.hdf5.0',fields)
                stars_forestdata = vpt.ReadForest(stf_folder + 'stars_forest.hdf5.0',fields)
            else:
                dm_forestdata = None
                stars_forestdata = None
            # Compute catalogue quantities for each ozy file (using MPI if possible)
            # TODO : Just send the snapshots required, not the full forest
            dm_forestdata = comm.bcast(dm_forestdata,root=0)
            stars_forestdata = comm.bcast(stars_forestdata,root=0)
            for snap in rank_snaps:
                # Get the index in the forest for the data (it's done for DM
                # but it should be the same for the stars!)
                ind = dm_forestdata[5].index(snap.snapname)
                snap.build_STF(skipran, dm_forestdata, stars_forestdata, ind, **kwargs)
    else:
        print('WARNING: The halo finder %s is not supported, pleae check!'%halofinder)
        exit

if __name__ == '__main__':
    print_art()
