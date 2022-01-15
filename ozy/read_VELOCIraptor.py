import sys
import os
import glob
import time
from tokenize import group
import numpy as np
import h5py
import copy

from ozy.group import create_new_group, grouptypes
from ozy.utils import remove_out_zoom

# Load VELOCIraptor python routines (written by the developers of VELOCIraptor)
# Load the cythonized code if compiled
try:
    import velociraptor_python_tools_cython as vpt
    print('Using cython VR+TF toolkit')
except:
    import velociraptor_python_tools as vpt
    print('Using python VR+TF toolkit')

plist_dict = dict(star='slist', bh='bhlist', dust='dlist', dm='dmlist', dm2='dm2list')


# Define properties of interest to be extracted from VELOCIraptor catalogues
# ID : Halo ID. ID = index of halo + 1 + TEMPORALHALOIDVAL * Snapshot_value, 
# giving a temporally unique halo id that can be quickly parsed for an index 
# and a snapshot number.
# hostHaloID : ID of the host field halo. If an object is a field halo, this is -1.
# Structuretype : Structure types contain information on how the object was found 
# and at what level in the subhalo hierarchy. Field halos are 10. Substructures 
# identified using the local velocity field are type 10+10=20, substructures 
# identified using cores are type 10+5=15. For structures found at level 2 
# (ie: subhalos within subhalos), the type offset is 20, and so on.
# numSubStruct : Number of substructures. Subhalos can have subsubhalos.
requestedfields=[
    'ID', 'hostHaloID',
    'Structuretype', 'Mass_tot', 'npart', 'numSubStruct'
    ]
# Mass_FOF : Total mass of particles in the FOF, Mfof. Is zero for substructure.
# Mass_200mean : Overdensity mass defined by mean matter density.
# Mass_200crit : Overdensity mass defined by critical density.
# R_size : Maximum distance of particles belonging exclusively to the object and the object’s centre.
# R_HalfMass : Half mass radius based on the Mass_tot.
# R_200mean : Same as above but for Mass_200mean.
# R_200crit : Same as above but for Mass_200crit.
# Xc,Yc,Zc : Coordinates of the center of mass.
# Xcminpot,Ycminpot,Zcminpot : Coordinates of the potential minimum.
# Vxc,Vyc,Vzc : Velocity of the center of mass.
# lamdda_B : Bullock et al. (2001) like spin parameter.
# Lx,Ly,Lz : Angular momentum vector about the COM using particles belonging to halo.
# RVmax_Lx,RVmax_Ly,RVmax_Lz : Same as above but at R(Vmax).
# sigV : Velocity dispersion based on the velocity dispersion tensor.
# RVmax_sigV : Same as above but at R(Vmax).
# Rmax : Radius of maximum circular velocity.
# Vmax : Maximum circular velocity based on particles belonging to the halo.
# cNFW : Concentration assuming NFW profile (Navarro, Frenk & Whitr 1997).
# Efrac : The fraction of particles formally bound.
# Epot : Total gravitational energy.
# Ekin : Total kinetic energy.
additionalrequestedfields=[
    'Mass_FOF', 'Mass_200mean', 'Mass_200crit',
    'R_size', 'R_HalfMass', 'R_200mean', 'R_200crit',
    'Xc', 'Yc', 'Zc',
    'Xcminpot', 'Ycminpot', 'Zcminpot',
    'VXc', 'VYc', 'VZc',
    'lambda_B',
    'Lx','Ly','Lz',
    'RVmax_Lx','RVmax_Ly','RVmax_Lz',
    'sigV', 'RVmax_sigV',
    'Rmax', 'Vmax',
    'cNFW',
    'Efrac', 'Epot', 'Ekin'
    ]

def make_walkable_tree(stf_folder,grouptype):
    """
    This function loads the raw tree files from TreeFrog and builds the head/tail, 
    root head/root tail info writing a unified walkable tree.
    """

    # Determine required halos file.
    if grouptype == 'halo':
        basename = 'dm'
        outputfname = stf_folder+'dm_walkable_tree'
    elif grouptype == 'galaxy':
        basename = 'stars'
        outputfname = stf_folder+'stars_walkable_tree'
    else:
        return
    
    # Get the names of the available snapshots
    snaps = glob.glob(stf_folder+basename+'.output_*.properties')
    snaps = [snap.split('.')[1] for snap in snaps]
    snaps = sorted(snaps, key=lambda x: int(x[-5:]))

    # Review whether requested structure catalogue folder exists
    if not os.path.exists(stf_folder):
        return

    # Requested halo fields
    requestedfields=[
        'ID', 'hostHaloID',
        ]

    # Define different types of input
    ASCIIINPUT = 0
    HDFINPUT = 2

    # TODO: Find a way to actually automatise this
    # Get the version of TreeFrog
    TFVERSION = 1.21
    HFNAME = 'VELOCIraptor'
    HFVERSION = 1.50

    # TODO: Could alter to have user indicate input type but currently assume all is HDF
    RAWTREEFORMAT = HDFINPUT
    RAWPROPFORMAT = HDFINPUT

    # Load the tree information stored in the file such as temporal halo ID,
    # number of snapshots search...
    hfile = h5py.File(stf_folder+basename+'.snapshot_000.VELOCIraptor.tree')
    numsnaps = hfile.attrs['Number_of_snapshots']
    TEMPORALHALOIDVAL = hfile.attrs['Temporal_halo_id_value']
    NSNAPSEARCH = hfile.attrs['Nsteps_search_new_links']
    TREEDIRECTION = hfile.attrs['Tree_direction']
    hfile.close()

    # Check that the tree files are updated
    if numsnaps != len(snaps):
        print('The present tree files are outdated compared to halo catalogues. Check!')
        exit

    # Number of particles used in the halo catalog
    # TODO: To be taken from configuration files
    NPARTTHRESHOLD = 100

    rawtreedata = None
    # Read raw descendant tree along with merits, don't reverse snap order
    snaptreelist = open(stf_folder+basename+'.snaptreelist.txt','w')
    for i in range(numsnaps):
        snaptreelist.write(stf_folder+basename+'.snapshot_%03d.VELOCIraptor\n'%i)
    snaptreelist.close()
    fname = stf_folder + basename + '.snaptreelist.txt'

    if (TREEDIRECTION == 1):
        rawtreedata = vpt.ReadHaloMergerTreeDescendant(fname, False, RAWTREEFORMAT, 1, True)
    elif (TREEDIRECTION == 0):
        print('Warning, progenitor based trees are less useful when building halo merger trees for SAMs.')
        rawtreedata = vpt.ReadHaloMergerTree(fname, RAWTREEFORMAT, 1, True)
    else:
        print('Full graphs to walkable trees not implemented yet.')
    
    print('Finished reading raw tree')
    numhalos = np.zeros(numsnaps,dtype=np.uint64)
    halodata = [dict() for i in range(numsnaps)]
    atime = np.zeros(numsnaps)
    for i in range(numsnaps):
        halodata[i],numhalos[i]=vpt.ReadPropertyFile(stf_folder+basename+'.'+snaps[i], 2, 0, 1, requestedfields)
        atime[i]=halodata[i]['SimulationInfo']['ScaleFactor']
        for key in halodata[i].keys():
            if (key == 'SimulationInfo' or key == 'UnitInfo' or key == "ConfigurationInfo"): continue
            if (halodata[i][key].dtype==np.float64):
                halodata[i][key] = np.array(halodata[i][key],dtype=np.float32)
    print('Finished getting halo IDs')

    # Produce head tail in ascending order
    start = time.process_time()
    print("Building head/tail ")
    vpt.BuildTemporalHeadTailDescendant(numsnaps, rawtreedata, numhalos, halodata,
        TEMPORALHALOIDVAL)
    print("Finished head/tail ", time.process_time()-start)

    #store the description
    DescriptionInfo={
            'Title':'Walkable Tree',
            'TreeBuilder' : {
                'Name' : 'TreeFrog',
                'Version' : TFVERSION,
                'Temporal_linking_length' : NSNAPSEARCH,
                'Temporal_halo_id_value' : TEMPORALHALOIDVAL,
                'Tree_direction' : TREEDIRECTION,
            },
            'HaloFinder' : {
                'Name' : HFNAME, 'Version' : HFVERSION,
                'Particle_num_threshold' : NPARTTHRESHOLD,
                },
            }
    vpt.WriteWalkableHDFTree(outputfname, snaps, rawtreedata, numhalos, halodata,
                                atime, DescriptionInfo)

def make_forest(stf_folder, grouptype, zoom):
    """
    This function is called after make_walkable_tree, as it loads the resulting walkable tree,
    associated halo properties and builds the sublinks, progenitor links and forest IDs.
    All is finally saved in a forest file which can be used by Ozymandias.
    """

    # Determine required halos file.
    if grouptype == 'halo':
        basename = 'dm'
        igas = 0
        istar = 0
        ibh = 0
        idm = 1
        walkabletreefile = stf_folder+'dm_walkable_tree'
        outputfname = stf_folder+'dm_forest'
    elif grouptype == 'galaxy':
        basename = 'stars'
        igas = 0
        istar = 1
        ibh = 0
        idm = 0
        walkabletreefile = stf_folder+'stars_walkable_tree'
        outputfname = stf_folder+'stars_forest'
    else:
        return
    
    # Get the names of the available snapshots
    snaps = glob.glob(stf_folder+basename+'.output_*.properties')
    snaps = [snap.split('.')[1] for snap in snaps]
    snaps = sorted(snaps, key=lambda x: int(x[-5:]))

    # Review whether requested structure catalogue folder exists
    if not os.path.exists(stf_folder):
        print('The given STF folder does not exist!')
        return
    # Check for walkable tree
    if not os.path.exists(walkabletreefile):
        print('The given walkable tree does not exist!')
        return

    # Define different types of input
    ASCIIINPUT = 0
    HDFINPUT = 2

    # TODO: Find a way to actually automatise this
    # Get the version of TreeFrog
    TFVERSION = 1.21
    HFNAME = 'VELOCIraptor'
    HFVERSION = 1.50

    # TODO: Could alter to have user indicate input type but currently assume all is HDF
    RAWTREEFORMAT=HDFINPUT
    RAWPROPFORMAT=HDFINPUT

    # Load the tree information stored in the walkable tree
    treedata,numsnaps = vpt.ReadWalkableHDFTree(walkabletreefile)
    numsnaps = treedata['Header']['Number_of_snapshots']
    TEMPORALHALOIDVAL = treedata['Header']['TreeBuilder']['Temporal_halo_id_value']
    NSNAPSEARCH = treedata['Header']['TreeBuilder']['Temporal_linking_length']
    TREEDIRECTION = treedata['Header']['TreeBuilder']['Tree_direction']
    NPARTTHRESHOLD = treedata['Header']['HaloFinder']['Particle_num_threshold']
    # Alias treedata to halodata as forest file will combine the data
    halodata = treedata
    # Store number of halos, scalefactors
    numhalos = np.zeros(numsnaps, dtype=np.int64)
    scalefactors = np.zeros(numsnaps)

    # Load halo properties file
    print('Loading halo properties ...')
    time1 = time.time()
    for i in range(numsnaps):
        fname = stf_folder+basename+'.'+snaps[i]
        halos, numhalos[i] = vpt.ReadPropertyFile(fname,RAWPROPFORMAT,0,0,requestedfields)
        scalefactors[i]=halos['SimulationInfo']['ScaleFactor']
        halodata[i].update(halos)
    print('Done', time.time() - time1)

    # Given walkable tree, determine the largest difference in snapshots
    # between an object and its head
    maxnsnapsearch=0
    for i in range(numsnaps):
        if (numhalos[i] == 0): continue
        headsnap = np.int64(halodata[i]['Head']/TEMPORALHALOIDVAL)
        maxs = np.max(headsnap-i)
        maxnsnapsearch = max(maxnsnapsearch, maxs)
    print('Walkable tree has maximum snaps search of ', maxnsnapsearch)
    sys.stdout.flush()

    # Generate subhalo links
    vpt.GenerateSubhaloLinks(numsnaps,numhalos,halodata)
    # Generate progenitor links
    vpt.GenerateProgenitorLinks(numsnaps,numhalos,halodata)

    # Building forest
    ireverseorder = False
    iverbose = 1
    iforestcheck = True # This uses extra compute and is generally unnecessary
    forestdata = vpt.GenerateForest(snaps, numhalos, halodata, scalefactors,
                                    maxnsnapsearch, ireverseorder, TEMPORALHALOIDVAL, 
                                    iverbose, iforestcheck)
    
    # Strip out simulation and unit data
    SimulationInfo = copy.deepcopy(halodata[0]['SimulationInfo'])
    UnitInfo = copy.deepcopy(halodata[0]['UnitInfo'])
    HaloFinderConfigurationInfo = copy.deepcopy(halodata[0]['ConfigurationInfo'])
    if SimulationInfo['Cosmological_Sim']:
        if not UnitInfo['Comoving_or_Physical']:
            # Convert period to comoving little h
            SimulationInfo['Period'] *= SimulationInfo['h_val']/SimulationInfo['ScaleFactor']
        del SimulationInfo['ScaleFactor']
        # Lets update the names for genesis style
        UnitInfo['Comoving_unit_flag'] = UnitInfo['Comoving_or_Physical']
    
    # Remove the unit info from each snapshot as it is unnecessary
    for i in range(numsnaps):
        del halodata[i]['SimulationInfo']
        del halodata[i]['UnitInfo']
        del halodata[i]['ConfigurationInfo']
    
    # Write the unified file that contains forest ids, properties,
    # description will have to be updated so as to use appropriate version numbers
    DescriptionInfo={
            'Title':'Forest',
            'TreeBuilder' : copy.deepcopy(treedata['Header']['TreeBuilder']),
            'HaloFinder' : copy.deepcopy(treedata['Header']['HaloFinder']),
            'Flag_subhalo_links':True, 'Flag_progenitor_links':True, 'Flag_forest_ids':True, 'Flag_sorted_forest':False,
            'ParticleInfo':{
                'Flag_dm':(idm==1), 'Flag_gas':(igas==1), 'Flag_star':(istar==1), 'Flag_bh':(ibh==1), 'Flag_zoom': zoom,
                'Particle_mass' : {'dm':-1, 'gas':-1, 'star':-1, 'bh':-1, 'lowres':-1}
                }
            }

    DescriptionInfo['HaloFinder']['Subhalo_Particle_num_threshold'] = NPARTTHRESHOLD
    DescriptionInfo['TreeBuilder']['Temporally_Unique_Halo_ID_Description'] = 'Snap_num*Temporal_linking_length+Index+1'
    
    vpt.WriteForest(outputfname, snaps, numhalos, halodata, forestdata, scalefactors,
                    DescriptionInfo, SimulationInfo, UnitInfo, HaloFinderConfigurationInfo)

    halolist = [None for i in range(numsnaps)]
    for i in range(numsnaps):
        halolist[i]=stf_folder+basename+'.'+snaps[i]
    vpt.ForestFileAddHaloData(outputfname, halolist, snaps, additionalrequestedfields)
    print('Forest file done.')

def read_forest_portion(obj, forestdata, ind, grouptype):
    """
    This function reads the details of the halo structures (and tree) found by VELOCIraptor
    and TreeFrog, and constructs the Group objects with the global quantities of interest.
    """

    halodata, nhalos, atime, simdata, unitdata, snapnames = forestdata
    mydata = halodata[ind]

    nobjs = 0
    nonzoom_halos = 0
    # Looping over halos...
    for i in range(0, nhalos):
        new_group = create_new_group(obj, grouptype)
        new_group.npart = int(mydata['npart'][i])
        new_group.nsub = int(mydata['numSubStruct'][i])
        # TODO: Get particle IDs from files
        if grouptype == 'halo':
            new_group.ndm = new_group.npart
            new_group.dmlist = []
        elif grouptype == 'galaxy':
            new_group.nstar = new_group.npart
            new_group.slist = []
        # Halo integers - structure hierarchy and tree structure
        new_group.ID = int(mydata['ID'][i])
        new_group.Descendant = int(mydata['Descendant'][i])
        new_group.DescendantSnap = int(mydata['DescendantSnap'][i])
        new_group.FinalDescendant = int(mydata['FinalDescendant'][i])
        new_group.FinalDescendantSnap = int(mydata['FinalDescendantSnap'][i])
        new_group.FirstProgenitor = int(mydata['FirstProgenitor'][i])
        new_group.FirstProgenitorSnap = int(mydata['FirstProgenitorSnap'][i])
        new_group.ForestID = int(mydata['ForestID'][i])
        new_group.ForestLevel = int(mydata['ForestLevel'][i])
        new_group.Head = int(mydata['Head'][i])
        new_group.HeadIndex = int(mydata['HeadIndex'][i])
        new_group.HeadSnap = int(mydata['HeadSnap'][i])
        new_group.NextProgenitor = int(mydata['NextProgenitor'][i])
        new_group.NextSubhalo = int(mydata['NextSubhalo'][i])
        new_group.NumProgen = int(mydata['Num_progen'][i])
        new_group.PreviousProgenitor = int(mydata['PreviousProgenitor'][i])
        new_group.PreviousSubhalo = int(mydata['PreviousSubhalo'][i])
        new_group.Progenitor = int(mydata['Progenitor'][i])
        new_group.ProgenitorSnap = int(mydata['ProgenitorSnap'][i])
        new_group.RightTail = int(mydata['RightTail'][i])
        new_group.RootHead = int(mydata['RootHead'][i])
        new_group.RootHeadIndex = int(mydata['RootHeadIndex'][i])
        new_group.RootHeadSnap = int(mydata['RootHeadSnap'][i])
        new_group.RootTail = int(mydata['RootTail'][i])
        new_group.RootTailIndex = int(mydata['RootTailIndex'][i])
        new_group.RootTailSnap = int(mydata['RootTailSnap'][i])
        new_group.Structuretype = int(mydata['Structuretype'][i])
        new_group.Tail = int(mydata['Tail'][i])
        new_group.TailIndex = int(mydata['TailIndex'][i])
        new_group.TailSnap = int(mydata['TailSnap'][i])
        new_group.HostHaloID = int(mydata['hostHaloID'][i])

        # Halo masses
        m = (mydata['Mass_tot'][i] * float(unitdata['Mass_unit_to_solarmass']))/obj.quantity(1.0, 'code_mass').in_units('Msun').d
        new_group.mass['total_stf'] = obj.quantity(m, 'code_mass')
        m = (mydata['Mass_200crit'][i] * float(unitdata['Mass_unit_to_solarmass']))/obj.quantity(1.0, 'code_mass').in_units('Msun').d
        new_group.mass['200crit'] = obj.quantity(m, 'code_mass')
        m = (mydata['Mass_200mean'][i] * float(unitdata['Mass_unit_to_solarmass']))/obj.quantity(1.0, 'code_mass').in_units('Msun').d
        new_group.mass['200mean'] = obj.quantity(m, 'code_mass')
        m = (mydata['Mass_FOF'][i] * float(unitdata['Mass_unit_to_solarmass']))/obj.quantity(1.0, 'code_mass').in_units('Msun').d
        new_group.mass['FOF'] = obj.quantity(m, 'code_mass')

        # Halo position
        Xc = (mydata['Xc'][i] * float(unitdata['Length_unit_to_kpc']))/obj.quantity(1.0, 'code_length').in_units('kpc').d
        Yc = (mydata['Yc'][i] * float(unitdata['Length_unit_to_kpc']))/obj.quantity(1.0, 'code_length').in_units('kpc').d
        Zc = (mydata['Zc'][i] * float(unitdata['Length_unit_to_kpc']))/obj.quantity(1.0, 'code_length').in_units('kpc').d
        new_group.position['COM'] = obj.array(np.array([Xc,Yc,Zc]), 'code_length')

        Xc = (mydata['Xcminpot'][i] * float(unitdata['Length_unit_to_kpc']))/obj.quantity(1.0, 'code_length').in_units('kpc').d
        Yc = (mydata['Ycminpot'][i] * float(unitdata['Length_unit_to_kpc']))/obj.quantity(1.0, 'code_length').in_units('kpc').d
        Zc = (mydata['Zcminpot'][i] * float(unitdata['Length_unit_to_kpc']))/obj.quantity(1.0, 'code_length').in_units('kpc').d
        new_group.position['minpot'] = obj.array(np.array([Xc,Yc,Zc]), 'code_length')

        # Halo velocity
        VXc = (mydata['VXc'][i] * float(unitdata['Velocity_unit_to_kms']))
        VYc = (mydata['VYc'][i] * float(unitdata['Velocity_unit_to_kms']))
        VZc = (mydata['VZc'][i] * float(unitdata['Velocity_unit_to_kms']))
        new_group.velocity['COM'] = obj.array(np.array([VXc,VYc,VZc]), 'km/s')
        Vmax = (mydata['Vmax'][i] * float(unitdata['Velocity_unit_to_kms']))
        new_group.velocity['Vmax'] = obj.quantity(Vmax, 'km/s')
        sigV = (mydata['sigV'][i] * float(unitdata['Velocity_unit_to_kms']))
        new_group.velocity['sigV'] = obj.quantity(sigV, 'km/s')
        Rmax_sigV = (mydata['Rmax_sigV'][i] * float(unitdata['Velocity_unit_to_kms']))
        new_group.velocity['Rmax_sigV'] = obj.quantity(Rmax_sigV, 'km/s')

        # Halo angular momenta
        Lx = mydata['Lx'][i] * float(unitdata['Velocity_unit_to_kms']) * float(unitdata['Mass_unit_to_solarmass']) * float(unitdata['Length_unit_to_kpc'])
        Ly = mydata['Ly'][i] * float(unitdata['Velocity_unit_to_kms']) * float(unitdata['Mass_unit_to_solarmass']) * float(unitdata['Length_unit_to_kpc'])
        Lz = mydata['Lz'][i] * float(unitdata['Velocity_unit_to_kms']) * float(unitdata['Mass_unit_to_solarmass']) * float(unitdata['Length_unit_to_kpc'])
        new_group.angular_mom['COM'] = obj.array(np.array([Lx,Ly,Lz]), 'Msun*kpc*km/s')

        Lx = mydata['RVmax_Lx'][i] * float(unitdata['Velocity_unit_to_kms']) * float(unitdata['Mass_unit_to_solarmass']) * float(unitdata['Length_unit_to_kpc'])
        Ly = mydata['RVmax_Ly'][i] * float(unitdata['Velocity_unit_to_kms']) * float(unitdata['Mass_unit_to_solarmass']) * float(unitdata['Length_unit_to_kpc'])
        Lz = mydata['RVmax_Lz'][i] * float(unitdata['Velocity_unit_to_kms']) * float(unitdata['Mass_unit_to_solarmass']) * float(unitdata['Length_unit_to_kpc'])
        new_group.angular_mom['RVmax'] = obj.array(np.array([Lx,Ly,Lz]), 'Msun*kpc*km/s')

        # Halo radii
        r = (mydata['R_200crit'][i] * float(unitdata['Length_unit_to_kpc']))/obj.quantity(1.0, 'code_length').in_units('kpc').d
        new_group.radius['200crit'] = obj.quantity(r, 'code_length')
        r = (mydata['R_200mean'][i] * float(unitdata['Length_unit_to_kpc']))/obj.quantity(1.0, 'code_length').in_units('kpc').d
        new_group.radius['200mean'] = obj.quantity(r, 'code_length')
        r = (mydata['R_HalfMass'][i] * float(unitdata['Length_unit_to_kpc']))/obj.quantity(1.0, 'code_length').in_units('kpc').d
        new_group.radius['HalfMass'] = obj.quantity(r, 'code_length')
        r = (mydata['R_size'][i] * float(unitdata['Length_unit_to_kpc']))/obj.quantity(1.0, 'code_length').in_units('kpc').d
        new_group.radius['size'] = obj.quantity(r, 'code_length')
        r = (mydata['Rmax'][i] * float(unitdata['Length_unit_to_kpc']))/obj.quantity(1.0, 'code_length').in_units('kpc').d
        new_group.radius['max'] = obj.quantity(r, 'code_length')
        r = (mydata['Rmax'][i] * float(unitdata['Length_unit_to_kpc']))/obj.quantity(1.0, 'code_length').in_units('kpc').d
        new_group.radius['max'] = obj.quantity(r, 'code_length')

        # Halo energies
        E = mydata['Ekin'][i] * float(unitdata['Mass_unit_to_solarmass']) * (float(unitdata['Velocity_unit_to_kms'])**2.0)
        new_group.energies['kinetic'] = obj.quantity(E, 'Msun * km**2 * s**-2')
        E = mydata['Epot'][i] * float(unitdata['Mass_unit_to_solarmass']) * (float(unitdata['Velocity_unit_to_kms'])**2.0)
        new_group.energies['potential'] = obj.quantity(E, 'Msun * km**2 * s**-2')
        new_group.BoundFrac = float(mydata['Efrac'][i])

        # Other properties
        new_group.cNFW = mydata['cNFW'][i]
        new_group.lambda_B = mydata['lambda_B'][i]

        # If the simulation is a zoom, check if we want to throw away
        # objects outside of the zoom regoion.
        add_group = True
        if obj.simulation.zoom:
            add_group = remove_out_zoom(obj, new_group)
        if add_group:
            if new_group._valid:
                obj.__dict__[grouptypes[grouptype]].append(new_group)
                nobjs += 1
        else:
            nonzoom_halos += 1

        if obj.simulation.zoom:
            print("Halos out of zoom region: "+str(nonzoom_halos))
        if grouptype == 'halo':
            obj.nhalos = nobjs
            print("Number of selected DM halos: "+str(obj.nhalos))
        elif grouptype == 'galaxy':
            obj.ngalaxies = nobjs
            print("Number of selected galaxies: "+str(obj.ngalaxies))