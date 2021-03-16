import os
import yt
import ozy
from ozy.progen import progen_build

class Snapshot(object):
    """Class for tracking paths and data for simulation snapshots.
    """
    def __init__(self, snapdir, snapname, snapindex):
        self.snapdir   = snapdir
        self.snapname  = snapname
        self.snapindex = snapindex
        self.snap      = '%s/%s%05d/info_%05d.txt' % (snapdir, snapname, snapindex, snapindex) # Assuming RAMSES convention of integer of length 5 for the snapindex
    
    def set_output_information(self, ds, prefix='ozy_', extension='hdf5'):
        
        self.outdir = '%s/Groups' % ('/' + os.path.join(*ds.fullpath.split('/')[0:-1]))
        
        self.outfile = '%s/%s%s%05d.%s' % (self.outdir, prefix, self.snapname.replace(self.snapname,''), self.snapindex, extension)
    
    def _make_output_dir(self):
        if not os.path.isdir(self.outdir):
            try:
                os.makedirs(self.outdir)
            except:
                pass
    def build_HaloMaker(self, skipran, **kwargs):
        
        if not os.path.isfile(self.snap):
            return
        ds = yt.load(self.snap)
        self.set_output_information(ds)
        
        if os.path.isfile(self.outfile) and skipran:
            return
        self._make_output_dir()
        
        obj = ozy.OZY(ds)
        obj.build_HaloMaker(**kwargs)
        obj.save(self.outfile)
        
        obj = None
        ds = None
        
def print_art():
    from ozy.__version__ import VERSION
    copywrite = '   (C) 2021 F. Rodriguez Montero'
    version   = '   Version %s' % VERSION
    
    art = """    )                                                       
    ( /(                                   (                   
    )\())      (        )       )          )\ )  (      )      
    ((_)\   (   )\ )    (     ( /(   (     (()/(  )\  ( /(  (   
    ((_)  )\ (()/(    )\  ' )(_))  )\ )   ((_))((_) )(_)) )\  
    / _ \ ((_) )(_)) _((_)) ((_)_  _(_/(   _| |  (_)((_)_ ((_) 
    | (_) ||_ /| || || '  \()/ _` || ' \))/ _` |  | |/ _` |(_-< 
    \___/ /__| \_, ||_|_|_| \__,_||_||_| \__,_|  |_|\__,_|/__/ 
                |__/                                            
                
    """
    print('\n%s\n%s\n%s\n' % (art, copywrite, version))

def drive(snapdirs, snapname, snapindexes, progen=False, skipran=False,
          build_HaloMaker=True, extension='hdf5', ozy_prefix='ozy_', **kwargs):
    """Driver function for running `ÒZYMANDIAS``on multiple snapshots.
    
    """
    # Just in the case we are giv a single simulation directory or
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
    snaps = []
    for snapdir in snapdirs:
        for snapindex in snapindexes:
            snaps.append(Snapshot(snapdir, snapname, snapindex))
    
    if build_HaloMaker:
        rank_snaps = snaps[rank::nprocs]
        for snap in rank_snaps:
            snap.build_HaloMaker(skipran, **kwargs)
    
    if progen:
        ozy.progen.run_progen(snapdirs, snapname, snapindexes, prefix=ozy_prefix, extension=extension, **kwargs)

if __name__ == '__main__':
    print_art()
