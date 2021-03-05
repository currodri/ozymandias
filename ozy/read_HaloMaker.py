import numpy as np
import os
from ozy.group import create_new_group, grouptypes
from ozy.utils import remove_out_zoom

plist_dict = dict( gas='glist', star='slist', bh='bhlist', dust='dlist', dm='dmlist', dm2='dm2list')

def read_HM(obj, grouptype):
    """Read found structures from HaloMaker files.
    
    This function reads the details of the halo structures found by the HaloMaker code
    (Tweed et al. 2009) and constructs the Group objects with the global quantities
    already computed.
    
    The structure of this function is based upon read_halos, written in Python by
    S. Martin-Alvarez.
    
    """
    # yT returns the fullpath something like '/mnt/extraspace/currodri/NUT/cosmoNUT/output_00013'
    # so need to get rid of the last bit to get the path to the simulation, not the snap.
    sim_folder = '/' + os.path.join(*obj.simulation.fullpath.split('/')[0:-1]) + '/'
    
    # Just to get the index of the snap so that we can find the specific brick file.
    snap_ID = obj.simulation.fullpath.split('/')[-1].split('_')[-1]
    
    # Determine required halos file.
    if grouptype == 'halo':
        haloM_folder = 'HaloMaker_DM/DMOnly/'
        read_file    = 'tree_bricks%03d' % int(snap_ID)
    elif grouptype == 'galaxy':
        haloM_folder = 'HaloMaker_stars/StarsOnly/'
        read_file    = 'tree_brick_starsub_%03d' % int(snap_ID)
    elif grouptype == 'cloud':
        haloM_folder = 'HaloMaker_stars/GasOnly/'
        read_file    = 'tree_brick_starsub_%03d' % int(snap_ID)
    else:
        return
    
    # Review whether requested structure catalogue exists.
    file_route = sim_folder + haloM_folder + read_file
    if not os.path.exists(file_route):
        return
    
    # If it does exist, find the cleaned version, if that option is used.
    clean_route = sim_folder + haloM_folder + 'clean_' + read_file
    if os.path.exists(clean_route) and obj.clean_brickfile:
        file_route = clean_route
    elif obj.clean_brickfile:
        clean_up_done = auto_cleanHM(sim_folder, haloM_folder)
        if clean_up_done:
            file_route = clean_route
    # Override if its a galaxy, since we need particles for progen search
    if grouptype == 'galaxy':
        file_route = sim_folder + haloM_folder + read_file
    # BEGIN - Open file to be read.
    HMfile = open(file_route, 'rb')
    
    # Reading header...
    nbodies = np.fromfile(file=HMfile, dtype=np.int32, count=3)
    nbodies = nbodies[1]
    realvals = np.fromfile(file=HMfile, dtype=np.float32, count=12)
    massp = realvals[1]
    aexp = realvals[4]
    omega_t = realvals[7]
    age_univ = realvals[10]
    nhalos = np.fromfile(file=HMfile, dtype=np.int32, count=4)
    nb_of_halos = nhalos[1]
    nb_of_subhalos = nhalos[2]
    
    if grouptype == 'halo':
        obj.nhalos = nb_of_halos + nb_of_subhalos
    elif grouptype == 'galaxy':
        obj.ngalaxies = nb_of_halos + nb_of_subhalos
    elif grouptype == 'cloud':
        obj.nclouds = nb_of_halos + nb_of_subhalos
    
    print( "Number of halos in output: " + str(nb_of_halos))
    print( "Number of subhalos in output: " + str(nb_of_subhalos))
    print( "Number of particles in output: " + str(nbodies))

    nobjs = 0
    nonzoom_halos = 0
    # Reading halos...
    for i in range(0, nb_of_halos + nb_of_subhalos):
        new_group = create_new_group(obj, grouptype)
        halo1 = np.fromfile(file=HMfile, dtype=np.int32, count=3)
        new_group.npart = int(halo1[1])
        if grouptype == 'halo':
            new_group.ndm = halo1[1]
        # TODO: Allow particle data to be store for each structure. 
        if not os.path.isfile(clean_route) or grouptype == 'galaxy':
            print('Number of star particles in galaxy: '+str(new_group.npart))
            ignore = np.fromfile(file=HMfile, dtype=np.int32, count=1)
            new_group.slist = np.fromfile(file=HMfile, dtype=np.int32, count=new_group.npart)
            ignore = np.fromfile(file=HMfile, dtype=np.int32, count=1)
        # Halo integers.
        tempR = np.fromfile(file=HMfile, dtype=np.int32, count=3)
        new_group.ID=tempR[1]
        tempR = np.fromfile(file=HMfile, dtype=np.int32, count=3)
        halo_tstep=tempR[1]
        tempR = np.fromfile(file=HMfile, dtype=np.int32, count=7)
        if grouptype == 'halo':
            new_group.level   = tempR[1]
            new_group.host    = tempR[2]
            new_group.hostsub = tempR[3]
            new_group.nsub    = tempR[4]
            new_group.nextsub = tempR[5]
        else:
            level   = tempR[1]
            host    = tempR[2]
            hostsub = tempR[3]
            nsub    = tempR[4]
            nextsub = tempR[5]
        # Halo total mass.
        tempR = np.fromfile(file=HMfile, dtype=np.float32, count=3)
        # Because of stupid HaloMaker, the units are in 10^11 Msun…
        m = (tempR[1] * 1e+11)/obj.yt_dataset.quan(1.0, 'code_mass').in_units('Msun').d
        new_group.mass['total'] = obj.yt_dataset.quan(m, 'code_mass')
        
        # Halo positions
        ignore = np.fromfile(file=HMfile, dtype=np.int32, count=1)
        tempR = np.fromfile(file=HMfile, dtype=np.float64, count=3)
        ignore = np.fromfile(file=HMfile, dtype=np.int32, count=1)
        new_group.position = obj.yt_dataset.arr(np.array([tempR[0]+0.5,tempR[1]+0.5,tempR[2]+0.5]), 'code_length')
        
        # Halo velocities
        ignore = np.fromfile(file=HMfile, dtype=np.int32, count=1)
        tempR = np.fromfile(file=HMfile, dtype=np.float64, count=3)
        ignore = np.fromfile(file=HMfile, dtype=np.int32, count=1)
        new_group.velocity = obj.yt_dataset.arr(np.array([tempR[0],tempR[1],tempR[2]]), 'code_velocity')
        
        # Halo angular momenta
        ignore = np.fromfile(file=HMfile, dtype=np.int32, count=1)
        tempR = np.fromfile(file=HMfile, dtype=np.float64, count=3)
        ignore = np.fromfile(file=HMfile, dtype=np.int32, count=1)
        new_group.angular_mom['total'] = obj.yt_dataset.arr(np.array([tempR[0],tempR[1],tempR[2]]), 
        'code_mass*code_length*code_velocity')
        
        # Halo radius + a,b,c
        tempR = np.fromfile(file=HMfile, dtype=np.float32, count=6)
        new_group.radius['total'] = obj.yt_dataset.quan(tempR[1], 'code_length')
        new_group.shape['total'] = np.array([tempR[2],tempR[3],tempR[4]])
        
        # Halo Kinetic,Potential,Thermal
        tempR = np.fromfile(file=HMfile, dtype=np.float32, count=5)
        new_group.energies['total_kinetic'] = obj.yt_dataset.quan(tempR[1], 'code_mass * code_velocity**2')
        new_group.energies['total_potential'] = obj.yt_dataset.quan(tempR[1], 'code_mass * code_velocity**2')
        
        
        # Halo Spin
        tempR = np.fromfile(file=HMfile, dtype=np.float32, count=3)
        if grouptype == 'halo':
            new_group.spin = tempR[1]
        
        # Halo virials [radius, mass, temperature, circular velocity]
        tempR = np.fromfile(file=HMfile, dtype=np.float32, count=6)
        new_group.virial_quantities['radius'] = obj.yt_dataset.quan(tempR[1], 'code_length')
        # Same here for the virial mass…
        m = (tempR[2] * 1e+11)/obj.yt_dataset.quan(1.0, 'code_mass').in_units('Msun').d
        new_group.virial_quantities['mass']   = obj.yt_dataset.quan(m, 'code_mass')
        new_group.virial_quantities['temperature']   = obj.yt_dataset.quan(tempR[3], 'code_temperature')
        new_group.virial_quantities['cvel']   = obj.yt_dataset.quan(tempR[4], 'code_velocity')
        
        # Halo profiles
        tempR = np.fromfile(file=HMfile, dtype=np.float32, count=4)
        if grouptype == 'halo':
            new_group.NFW_rho0 = tempR[1]
            new_group.NFW_r_c = tempR[2]
        # If the simulation is a zoom, check if we want to throw away
        # objects outside of the zoom regoion.
        add_group = True
        if obj.simulation.zoom:
            add_group = remove_out_zoom(obj, new_group)
        if add_group:
            obj.__dict__[grouptypes[grouptype]].append(new_group)
            nobjs += 1
        else:
            nonzoom_halos += 1
        
    halo1 = np.fromfile(file=HMfile, dtype=np.int32, count=3)
        
    HMfile.close()
    if obj.simulation.zoom:
        print("Halos out of zoom region: "+str(nonzoom_halos))

    if grouptype == 'halo':
        obj.nhalos = nobjs
    elif grouptype == 'galaxy':
        obj.ngalaxies = nobjs
    elif grouptype == 'cloud':
        obj.nclouds = nobjs
    
def auto_cleanHM(sim_folder, haloM_folder):
    """This function installs and automatically cleans the HaloMaker files for faster execution.
    """
    
    # TODO: The HalosExtractor should be a submodule of OZYMANDIAS, and compiled during
    # the installation of the package.

    clean_up_done = False
    HalosExtractorRoute = '/mnt/zfsusers/currodri/bin/' # This requires the "make" of the Fortran code first.
    if (os.path.isfile(HalosExtractorRoute+"HaloExtractor.out")):        
        print("I am trying to clean the halos in the received folder to accelerate")
        print("execution. This is done once per simulation and really worth not   ")
        print("having to read the particles")
        os.system(str(HalosExtractorRoute)+"./HaloExtractor.out "+str(sim_folder+haloM_folder))
        clean_up_done = True
    else:
        print("I have tried to clean your halo files for faster execution but")
        print("I could not find HalosExtractor.out. Consider installing this")
        print("for faster execution")
    return clean_up_done