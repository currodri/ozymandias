import numpy as np
import h5py
import os
import ozy
from ozy.group import grouptypes
from joblib import Parallel, delayed
from scipy import stats

def run_progen(snapdirs, snapname, snapindexes, prefix='ozy_', extension='hdf5', **kwargs):
    """Function to run progenitor/descendant finder in specified snapshot folders in a given directory.
    """
    from ozy.driver import Snapshot

    # Find existing snapshots in snapdirs
    if isinstance(snapdirs, str):
        snapdirs = [snapdirs]
    if isinstance(snapindexes, int):
        snapindexes = [int]
    
    # Make sure that snapindexes are sorted
    snapindexes.sort(reverse=True)
    snaps = []
    for snapdir in snapdirs:
        for snapindex in snapindexes:
            snaps.append(Snapshot(snapdir, snapname, snapindex))
    
    verified_snaps = []
    missing_snaps = []
    missing = ''
    for isnap, snap in enumerate(snaps):
        fullname = snap.snap
        if not os.path.isfile(fullname):
            missing_snaps.append(snap)
            missing = missing+' %05d'%(snapindexes[isnap])
            continue
        snap.outfile = ozy_filename(snap, prefix, extension)
        if not os.path.isfile(snap.outfile):
            missing_snaps.append(snap)
            missing = missing+' %05d'%(snapindexes[isnap])
            continue
        f = h5py.File(snap.outfile, 'r')
        if not '/halo_data' in f:
            missing_snaps.append(snap)
            missing = missing+' %05d'%(snapindexes[isnap])
            f.close()
            continue
        verified_snaps.append(snap)
        f.close()
    
    if len(missing_snaps) > 0:
        print('Missing output/ozymandias file, or no halo_data for: %s'%missing)
    
    # Collect pairs of snapshot names over which to run progen
    progen_pairs = []
    for i in range(0, len(verified_snaps)-1):
        progen_pairs.append((verified_snaps[i], verified_snaps[i+1]))
    
    # Loop over pairs, find progenitors
    for pair in progen_pairs:
        snap_current = pair[0]
        snap_progens = pair[1]

        if snap_current.snapindex < snap_progens.snapindex: 
            print('Progen: Finding descendants of snap %d in snap %d'%(snap_current.snapindex,snap_progens.snapindex))
        else:
            print('Progen: Finding progenitors of snap %d in snap %d'%(snap_current.snapindex,snap_progens.snapindex))
        
        obj_current = ozy.load(ozy_filename(snap_current, prefix, extension))
        obj_progens = ozy.load(ozy_filename(snap_progens, prefix, extension))

        progen_build(obj_current, obj_progens, ozy_filename(snap_current, prefix, extension), snap_dir = snapdirs[0], **kwargs)

def progen_build(obj_current, obj_target, ozy_file, snap_dir=None, data_type='galaxy', part_type='star', recompute=False, 
                    save=True, n_most=1, min_in_common=0.1, nproc=1):
    """Function to find the most massive progenitor of each Ozy group in obj_current in the previous snapshot.
    """

    if obj_current.simulation.redshift > obj_target.simulation.redshift:
        index_name = 'descend_'+data_type+'_'+part_type
    else:
        index_name = 'progen_'+data_type+'_'+part_type

    if not recompute and check_if_progen_present(ozy_file, index_name):
        print('%s data already present; returning data (set recompute=True to recompute)!' % (index_name))
        f = h5py.File(ozy_file,'r')
        prog_indexes = f['tree_data/%s'%index_name]
        return np.asarray(prog_indexes)    

    ng_current, pid_current, gid_current, pid_hash = collect_group_IDs(obj_current, data_type, part_type)
    ng_target, pid_target, gid_target, _ = collect_group_IDs(obj_target, data_type, part_type)

    if ng_current == 0 or ng_target == 0:
        print('No %s found in current ozy/target file (%d/%d) -- exiting progen_finder'%(data_type,ng_current,ng_target))
        return None
    # Now run the finder of progenitor indexes
    prog_indexes = find_progens(pid_current, pid_target, gid_current, gid_target, pid_hash, n_most=n_most, min_in_common=min_in_common, nproc=nproc)

    if save:
        write_progens(obj_current, np.array(prog_indexes).T, ozy_file, index_name, obj_target.simulation.redshift)
    
    return prog_indexes

def find_progens(pid_current, pid_target, gid_current, gid_target, pid_hash, n_most=1, min_in_common=0.1, nproc=1):
    """Distribute search of progenitors/descendants for each target galaxy.
    """
    # Sort the progenitors IDs and object numbers for faster searching
    isort_target = np.argsort(pid_target)
    pid_target = pid_target[isort_target]
    gid_target = gid_target[isort_target]
    ngroups_curr = len(pid_hash) - 1
    
    # Loop over current_objects to find progens for each
    if nproc > 1:
        prog_index_tmp = Parallel(n_jobs=nproc)(delayed(_find_target_group)(pid_current[pid_hash[ig]:pid_hash[ig+1]],pid_target,gid_target,min_in_common) for ig in range(ngroups_curr))
        prog_index_tmp = np.array(prog_index_tmp,dtype=int)
        prog_index = np.array(prog_index_tmp.T[0],dtype=int)
        prog_index2 = np.array(prog_index_tmp.T[1],dtype=int)
    else:
        prog_index = np.zeros(ngroups_curr,dtype=int)
        prog_index2 = np.zeros(ngroups_curr,dtype=int)
        for ig in range(ngroups_curr):
            prog_index[ig],prog_index2[ig] = _find_target_group(pid_current[pid_hash[ig]:pid_hash[ig+1]],pid_target,gid_target,min_in_common)
    
    if n_most == 1:
        return prog_index
    elif n_most == 2:
        return prog_index, prog_index2
    else:
        print('n_most=%d but must be 1 or 2; using 1'%n_most)
        return prog_index

def _find_target_group(pid_curr, pid_targ, gid_targ, min_in_common):
    """This is the main function which proceeds with the search of the most and second most particle ID's in common with pid_curr.
    """
    # Get where pid_curr would be in the order of pid_targ
    targ_ind = np.searchsorted(pid_targ, pid_curr)
    targ_ind = np.where(targ_ind==len(pid_targ),len(pid_targ)-1,targ_ind)
    ig_matched = np.where(pid_targ[targ_ind]==pid_curr,gid_targ[targ_ind],-1)
    ig_matched = ig_matched[ig_matched>=0]
    unique, counts = np.unique(ig_matched,return_counts=True)

    
    if len(ig_matched)>int(min_in_common*len(pid_curr)):
        # 1 Find target group index with most matches
        modestats = stats.mode(ig_matched)

        # 2 Store target group numbers
        prog_index_ig = modestats[0][0]

        # 3 Remove the first-most common group, recompute mode
        ig_matched = ig_matched[(ig_matched!=prog_index_ig)]
    else:
        # Set progen index to -1 when not found
        prog_index_ig = -1
    
    if len(ig_matched) > 0:
        # Now find it for the second-most matches
        modestats = stats.mode(ig_matched)
        prog_index_ig2 = modestats[0][0]
    else:
        prog_index_ig2 = -1
    
    return prog_index_ig, prog_index_ig2


def collect_group_IDs(obj, data_type, part_type):
    
    from ozy.read_HaloMaker import plist_dict

    if data_type == 'halo':
        part_list = 'obj.halos[i].%s'%plist_dict[part_type]
        ngroups = len(obj.halos)
    elif data_type == 'galaxy':
        part_list = 'obj.galaxies[i].%s'%plist_dict[part_type]
        ngroups = len(obj.galaxies)
    
    # Count number of total particles in groups
    npart = 0
    for i in range(ngroups):
        mylist = eval(part_list)
        npart += len(mylist)
    
    # Fill particle and group ID lists
    pids = np.zeros(npart,dtype=np.int64)
    gids = np.zeros(npart,dtype=np.int32)
    pid_hash = np.zeros(ngroups,dtype=np.int64)
    count = 0
    for i in range(0, ngroups):
        mylist = eval(part_list)
        pids[count:count+len(mylist)] = abs(mylist) # Star IDs are negative if they have not gone SN
        gids[count:count+len(mylist)] = np.full(len(mylist), i)
        pid_hash[i] = count
        count += len(mylist)
    pid_hash = np.append(pid_hash,npart+1)

    return ngroups, pids, gids, pid_hash


def write_progens(obj, data, ozy_file, index_name, redshift):
    """This function adds the results of the progen code to the original Ozy files.
    """

    f = h5py.File(ozy_file, 'r+')
    # Check if progen data is already present
    if check_if_progen_present(ozy_file, index_name):
        del f['tree_data/%s' % (index_name)]
        print('Overwriting %s info in tree_data' % (index_name))
    else:
        print('Writing %s info into tree_data'%(index_name))
    
    # Create group in HDF5 file
    try:
        tree = f.create_group('tree_data')
    except:
        tree = f['tree_data']
    progens = tree.create_dataset('%s' % (index_name), data=data, compression=1)
    tree.attrs[('z_'+index_name).encode('utf8')] = redshift
    f.close()
    return

def check_if_progen_present(ozy_file, index_name):
    """Check Ozy file for progren indexes.
    """

    f = h5py.File(ozy_file,'r')
    present = False
    if 'tree_data/%s' % (index_name) in f:
        present = True
    f.close()
    return present

def get_progen_redshift(ozy_file, index_name):
    """Return redshift of progenitors/descendants currently stored in tree_data.
    """

    f = h5py.File(ozy_file,'r')
    try:
        tree = f['tree_data']
    except:
        return -1
        print('Progen data %s does not exist in %s'%(index_name,ozy_file))
    z = tree.attrs['z_'+index_name]
    return z

def wipe_progen_info(ozy_file, index_name=None):
    """Get rid of any progen data in Ozy file.
    """
    f = h5py.File(ozy_file, 'r+')
    for dataset in f.keys():
        for name in f[dataset].keys():
            if index_name is None:
                if 'progen' in name or 'descend' in name:
                    print('Deleting %s from %s in %s'%(name,dataset,ozy_file))
                    del f[dataset][name]
            else:
                if index_name in name:
                    print('Deleting %s from %s in %s'%(name,dataset,ozy_file))
                    del f[dataset][name]
    f.close()

def ozy_filename(snap, prefix, extension):
    """Return full Ozy filename so that it can be loaded properly.
    """
    outdir = '%s/Groups' % (snap.snapdir)
    complete_path = '%s/%s%s%05d.%s' % (outdir, prefix, snap.snapname.replace(snap.snapname,''), snap.snapindex, extension)
    
    return complete_path