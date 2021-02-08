import numpy as np
import os
from ozy.group import create_new_group, grouptypes
from ozy.utils import remove_out_zoom

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
    sim_folder = '/' + os.path.join(*obj.simulation.fullpath.split('/')[0:-1])
    
    # Just to get the index of the snap so that we can find the specific brick file.
    snap_ID = obj.simulation.fullpath.split('/')[-1].split('_')[-1]
    
    # Determine required halos file.
    if grouptype == 'halo':
        haloM_folder = 'HaloMaker_DM/DMOnly/'
        read_file    = 'tree_bricks' + snap_ID
    elif grouptype == 'galaxy':
        haloM_folder = 'HaloMaker_stars/StarsOnly'
        read_file    = 'tree_starsub_' + snap_ID
    elif grouptype == 'cloud':
        haloM_folder = 'HaloMaker_stars/GasOnly'
        read_file    = 'tree_starsub_' + snap_ID
    else:
        return
    
    # Review whether requested structure catalogue exists.
    file_route = sim_folder + haloM_folder + read_file
    if not os.path.exists(file_route):
        return
    
    # If it does exist, find the cleaned version, if that option is used.
    clean_route = sim_folder + haloM_folder + 'clean_' + read_file
    if os.path.exists(file_route) and obj.clean_brickfile:
        file_route = clean_route
    elif obj.clean_brickfile:
        clean_up_done = auto_cleanHM(sim_folder, haloM_folder)
        if clean_up_done:
            file_route = clean_route
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
    
    # Reading halos...
    for i in range(0, nb_of_halos + nb_of_subhalos):
        new_group = create_new_group(obj, grouptype)
        halo1 = np.fromfile(file=HMfile, dtype=np.int32, count=3)
        new_group.npart = halo1[1]
        # TODO: Allow particle data to be store for each structure. 
        if not os.path.exists(file_route) and not obj.clean_brickfile:
            ignore = np.fromfile(file=HMfile, dtype=np.int32, count=1)
            for j in range(0, new_group.npart):
                partID = np.fromfile(file=HMfile, dtype=np.int32, count=1)
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
        # Halo total mass.
        tempR = np.fromfile(file=HMfile, dtype=np.float32, count=3)
        new_group.mass = tempR[1]
        
        # Halo positions
        ignore = np.fromfile(file=HMfile, dtype=np.int32, count=1)
        tempR = np.fromfile(file=HMfile, dtype=np.float64, count=3)
        ignore = np.fromfile(file=HMfile, dtype=np.int32, count=1)
        new_group.position = np.array([tempR[0]+0.5,tempR[1]+0.5,tempR[2]+0.5])
        
        # Halo velocities
        ignore = np.fromfile(file=HMfile, dtype=np.int32, count=1)
        tempR = np.fromfile(file=HMfile, dtype=np.float64, count=3)
        ignore = np.fromfile(file=HMfile, dtype=np.int32, count=1)
        new_group.velocity = np.array([tempR[0],tempR[1],tempR[2]])
        
        # Halo angular momenta
        ignore = np.fromfile(file=HMfile, dtype=np.int32, count=1)
        tempR = np.fromfile(file=HMfile, dtype=np.float64, count=3)
        ignore = np.fromfile(file=HMfile, dtype=np.int32, count=1)
        new_group.angular_mom = np.array([tempR[0],tempR[1],tempR[2]])
        
        # Halo radius + a,b,c
        tempR = np.fromfile(file=HMfile, dtype=np.float32, count=6)
        new_group.radius = tempR[1]
        new_group.shape = np.array([tempR[2],tempR[3],tempR[4]])
        
        # Halo Kinetic,Potential,Thermal
        tempR = np.fromfile(file=HMfile, dtype=np.float32, count=5)
        new_group.energies = np.array([tempR[1],tempR[2],tempR[3]])
        
        # Halo Spin
        tempR = np.fromfile(file=HMfile, dtype=np.float32, count=3)
        if grouptype == 'halo':
            new_group.spin = tempR[1]
        
        # Halo virials [radius, mass, temperature, circular velocity]
        tempR = np.fromfile(file=HMfile, dtype=np.float32, count=6)
        new_group.virial_radius = tempR[1]
        new_group.virial_mass = tempR[2]
        new_group.virial_temp = tempR[3]
        new_group.virial_cvel = tempR[4]
        
        # If the simulation is a zoom, check if we want to throw away
        # objects outside of the zoom regoion.
        add_group = True
        if obj.sim_attributes.zoom:
            add_group = remove_out_zoom(obj, new_group)
        if add_group:
            obj.__dict__[grouptypes[grouptype]].append(new_group)
        
        halo1 = np.fromfile(file=HMfile, dtype=np.int32, count=3)
        
    HMfile.close()
    
def auto_cleanHM(sim_folder, haloM_folder):
    """This function installs and automatically cleans the HaloMaker files for faster execution.
    """
    
    # TODO: The HalosExtractor should be a submodule of OZYMANDIAS, and compiled during
    # the installation of the package.
    
    HalosExtractorRoute = '~/bin' # This requires the "make" of the Fortran code first.
    if (os.path.isfile(HalosExtractorRoute+"HalosExtractor.out")):        
        print("I am trying to clean the halos in the received folder to accelerate")
        print("execution. This is done once per simulation and really worth not   ")
        print("having to read the particles")
        os.system(str(HalosExtractorRoute)+"./HalosExtractor.out "+str(sim_folder+haloM_folder))
    else:
        print("I have tried to clean your halo files for faster execution but")
        print("I could not find HalosExtractor.out. Consider installing this")
        print("for faster execution")