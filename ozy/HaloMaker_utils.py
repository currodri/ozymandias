import os
import numpy as np
from __init__ import OZYPATH
from hutils import py_halo_utils as phu
from ozy.group import create_new_group, grouptypes

class hmCatalogue(object):

    def __init__(self,Obj,RunDir,OutID,HaloDir='Halos',
                load=False):
        self.Obj = Obj # The OZY snapshot
        self.RunDir = RunDir # Where the simulation is located
        self.OutID = OutID # The output ID of the snapshot of interest
        self.OutDir = os.path.join(RunDir,'output_%5.5i'%self.OutID)
        self.info = Obj._info # The simulation info
        self.is_zoom = Obj.simulation.zoom # Whether the simulation is a zoom simulation
        self.HaloDir = os.path.join(RunDir,HaloDir,'output_%5.5i'%self.OutID) 
        self.TreeFile = os.path.join(self.HaloDir, 'tree_bricks%3.3i' % (self.OutID))
        self.is_galaxies = False # Root hmCatalogue is not for galaxies
        if 'tracer_particles' in Obj._kwargs:
            self.SimulationHasTracerParticles = Obj._kwargs['tracer_particles']
        else:
            self.SimulationHasTracerParticles = False
        if load:
            self.load_catalogue()

    def setup_HaloFinderRun(self, run=False, nvoisins=32, rhot=80., nhop=16, npart=20,
                                fudgepsilon=1e-5, alphap=1., method='MSM',cenmethod='com',
                                inputFile=None, ffile=None, runCmd=None, omega_b=False):
        """Set up the HaloFinder run with the specified parameters."""

        # 1. Check if the HaloMaker directory exists, if not create it
        if inputFile==None: inputFile = os.path.join(self.HaloDir, 'input_HaloMaker.dat')
        if ffile==None: ffile = os.path.join(self.HaloDir, 'inputfiles_HaloMaker.dat')
        if self.is_zoom:
            if runCmd==None: runCmd = OZYPATH + '/HaloFinder/program/HaloFinder_Zoom ./ > runHalos.log'
        else:
            if runCmd==None: runCmd = OZYPATH + '/HaloFinder/program/HaloFinder ./ > runHalos.log'
        if not os.path.exists(self.HaloDir):
            os.makedirs(self.HaloDir)
        
        # 2. Write the parameter file for HaloFinder
        boxlen_cMpc = self.info['boxlen'] * self.info['unit_l'] / 3.08e24 / self.info['aexp']
        f = open(inputFile, 'w')
        f.write("af              = 1.0        ! Today expansion factor == 1+zi \n")
        f.write("lbox            = %f         ! comoving size of the box (in Mpc)  \n"%(boxlen_cMpc))
        f.write("H_f             = %f         ! Hubble constant (in Km/s/Mpc) at z=0  \n"%(self.info['H0']))
        f.write("omega_f         = %f         ! Omega matter at z=0 (ex: 0.3)  \n"%(self.info['omega_m']))
        if omega_b and self.info['omega_b'] != 0.0:
            f.write("omega_b         = %f     ! Omega baryons at z=0 (ex: 0.3) \n"%(self.info['omega_b']))
        f.write("lambda_f        = %f         ! Omega lambda at z=0 (ex: 0.7)  \n"%(self.info['omega_l']))
        f.write("FlagPeriod      = 1          ! pericidicity flag choose 1 for periodic boundary conditions 0 else  \n")
        f.write("npart           = %i         ! Minimum number of particles per halo (ex:10)  \n"%(npart))
        f.write("method          = %s         ! choose between FOF, HOP, DPM, MHM (view ReadMe for details).  \n"%(method))
        if cenmethod == 'com':
            f.write("cdm             = .false.    ! use center of mass instead of densest particle.  \n")
            f.write("ssm             = .false.    ! use shrinking sphere method to find centre of halo.  \n")
        elif cenmethod == 'ssm':
            f.write("cdm             = .true.     ! use center of mass instead of densest particle.  \n")
            f.write("ssm             = .true.     ! use shrinking sphere method to find centre of halo.  \n")
        else:
            f.write("cdm             = .true.     ! use center of mass instead of densest particle.  \n")
            f.write("ssm             = .false.    ! use shrinking sphere method to find centre of halo.  \n")
        f.write("b               = 0.2        ! fof parameter b (usually b=0.2)  \n")
        f.write("nsteps          = 1          ! Number of time steps to analyse  \n")
        f.write("nvoisins        = %i         ! parameter for adaptahop (usually 20) \n"%(nvoisins))
        f.write("nhop            = %i         ! parameter for adaptahop (usually 20) \n"%nhop)
        f.write("rhot            = %f         ! parameter for adaptahop (80 coreespond to b = 0.2) \n"%(rhot))
        f.write("fudge           = 4.         ! parameter for adaptahop (usually 4.)  \n")
        f.write("fudgepsilon     = %f         ! parameter for adaptahop (usually 0.05) \n"%(fudgepsilon))
        f.write("alphap          = %f         ! parameter for adaptahop (usually 1.) \n"%(alphap))
        f.write("verbose         = .false.    ! verbose parameter for both halo finder \n")
        f.write("megaverbose     = .false.    ! parameter for adaptahop \n")
        if self.SimulationHasTracerParticles:
            f.write("SimulationHasTracerParticles = %r \n"%(self.SimulationHasTracerParticles))
        f.close()

        # 3. Write the input files for HaloFinder
        f = open(ffile, 'w')
        f.write("'%s/'  Ra3  1  %5.5i \n"%(self.OutDir,self.OutID))
        f.close()

        # 4. If requested, run HaloFinder
        if run and not os.path.exists(self.TreeFile):
            here = os.getcwd()
            os.chdir(self.HaloDir)
            os.system(runCmd)
            os.chdir(here)
    
    def load_catalogue(self):

        # 1. Check if the HaloMaker output directory exists already
        if not os.path.exists(self.HaloDir):
            raise FileNotFoundError(f"HalosFinder output directory {self.HaloDir} does not exist. Please run setup_HaloFinderRun first.")
        if not os.path.exists(self.TreeFile):
            raise FileNotFoundError(f"Tree file {self.TreeFile} does not exist. Please run setup_HaloFinderRun first.")
        
        # 2. Get the number of halos and subhalos
        nh,ns = phu.get_nb_halos(self.TreeFile)
        self.nhalos = nh
        self.nsubhalos = ns
        if nh == 0: return

        # 3. If halos are found, load the halo catalogue
        halos = np.zeros((nh+ns,32),dtype=np.float32,order='F')
        if self.is_galaxies  == False:
            grouptype = 'halo'
            if self.is_zoom:
                phu.read_all_halos_with_contam(self.TreeFile,halos,nh+ns)
            else:
                phu.read_all_halos(self.TreeFile,halos,nh+ns)
        else:
            grouptype = 'galaxy'
            phu.read_all_galaxies(self.TreeFile,halos,nh+ns)

        # 4. Create the groups and fill them with the halo data
        pos_offset = 0.5 * self.info['boxlen'] * self.info['unit_l'] / 3.08e24 # Half the boxlen in pMpc
        hm_correction_fact = 3.08e24 / self.Obj.quantity(1.0, 'Mpc').to('cm').d # To correct for the HM incorrect Mpc units
        nselected_halos, nselected_subhalos = 0, 0
        for i in range(0, nh+ns):
            new_group = create_new_group(self.Obj,grouptype)

            # Halo hierarchy
            new_group.npart = int(halos[i,0]) # Number of particles in the halo
            new_group.ID = int(halos[i,1]) # Halo ID
            new_group.tstep = int(halos[i,2]) # Time step of the snapshot (for TreeMaker)
            new_group._index = i # Index of the halo in the catalogue
            new_group.level = int(halos[i,3]) # Level of the halo in the hierarchy
            new_group.host = int(halos[i,4]) # Host halo ID (if any)
            new_group.hostsub = int(halos[i,5]) # Host subhalo ID (if any)
            new_group.nsub = int(halos[i,6])# Number of subhalos in this halo
            new_group.nextsub = int(halos[i,7]) # Next subhalo ID (if any)

            # Halo basic properties
            new_group.mass['total'] = self.Obj.quantity(halos[i,8]*1e11, 'Msun') # Halo mass (stupid HM with 10^11 Msun units...)
            new_group.position = self.Obj.array(np.array([halos[i,9]+pos_offset,
                                                           halos[i,10]+pos_offset,
                                                            halos[i,11]+pos_offset]), 'Mpc') # Halo position
            new_group.position = new_group.position * hm_correction_fact # Correct for HM's Mpc units
            new_group.velocity = self.Obj.array(np.array([halos[i,12], halos[i,13], halos[i,14]]), 'km/s') # Halo velocity
            new_group.angular_mom['total'] = self.Obj.array(np.array([halos[i,15]*1e11, halos[i,16]*1e11, halos[i,17]*1e11]), 'Msun*km/s*Mpc') # Halo angular momentum
            new_group.radius['total'] = self.Obj.quantity(halos[i,18] * hm_correction_fact, 'Mpc') # Halo radius
            new_group.shape['total'] = self.Obj.array(np.array([halos[i,19] * hm_correction_fact,
                                                                halos[i,20] * hm_correction_fact,
                                                                halos[i,21] * hm_correction_fact]), 'Mpc') # Halo shape (a,b,c)
            
            # Halo energies
            new_group.energies['total_kinetic'] = self.Obj.quantity(halos[i,22]*1e11, 'Msun*km**2/s**2') # Halo kinetic energy
            new_group.energies['total_potential'] = self.Obj.quantity(halos[i,23]*1e11, 'Msun*km**2/s**2') # Halo potential energy
            new_group.energies['total_binding'] = self.Obj.quantity(halos[i,24]*1e11, 'Msun*km**2/s**2') # Halo binding energy (kinetic + potential)
            new_group.spin = halos[i,25] # Halo spin parameter

            # Halo virial properties
            new_group.virial_quantities['radius'] = self.Obj.quantity(halos[i,26] * hm_correction_fact, 'Mpc') # Halo virial radius
            new_group.virial_quantities['mass'] = self.Obj.quantity(halos[i,27]*1e11, 'Msun') # Halo virial mass
            new_group.virial_quantities['temperature'] = self.Obj.quantity(halos[i,28], 'K') # Halo virial temperature
            new_group.virial_quantities['cvel'] = self.Obj.quantity(halos[i,29], 'km/s') # Halo virial circular velocity

            # Halo profile parameters (NFW or Isothermal Sphere)
            new_group.rho_0 = self.Obj.quantity(halos[i,30]*1e11/hm_correction_fact**3., 'Msun/Mpc**3') # Halo profile central density
            new_group.r_c = self.Obj.quantity(halos[i,31] * hm_correction_fact, 'Mpc') # Halo profile scale radius

            print(grouptype, new_group.ID, new_group.npart, new_group.mass['total']/1e6, new_group.virial_quantities['radius'].to('kpc'), new_group.radius['total'].to('kpc'))
            # Halo contamination (for zoom simulations)
            if self.is_zoom:
                new_group.contamination = halos[i,32]
            
            if new_group._valid:
                self.Obj.__dict__[grouptypes[grouptype]].append(new_group)
                if new_group.host == new_group.ID:
                    nselected_halos += 1
                else:
                    nselected_subhalos += 1

        # 5. Add up the number of halos and subhalos
        if self.is_galaxies:
            self.Obj.ngalaxies = nselected_halos + nselected_subhalos
            self.Obj.nsatellites = nselected_subhalos
            print(f"Loaded {self.Obj.ngalaxies} galaxies and {self.Obj.nsatellites} satellites from {self.TreeFile}.")
        else:
            self.Obj.nhalos = nselected_halos + nselected_subhalos
            self.Obj.nsubhalos = nselected_subhalos
            print(f"Loaded {self.Obj.nhalos} halos and {self.Obj.nsubhalos} subhalos from {self.TreeFile}.")
        
class galaxyCatalogue(hmCatalogue):
    """Class for galaxy catalogues, inheriting from hmCatalogue."""
    
    def __init__(self, Obj, RunDir, OutID, HaloDir='Galaxies', load=False):
        super().__init__(Obj, RunDir, OutID, HaloDir, load)
        self.is_galaxies = True # This is a galaxy catalogue
        self.TreeFile = os.path.join(self.HaloDir, 'tree_bricks_stars%3.3i' % (self.OutID))
        if load:
            self.load_catalogue()

    def setup_GalFinderRun(self, run=False, nvoisins=10, rhot=1000., nhop=10, npart=10,
                                fudgepsilon=1e-5, alphap=1., method='MSM'):
        """Set up the GalaxyFinder run with the specified parameters."""
        # Call the parent method with GalaxyFinder specific parameters
        inputFile = os.path.join(self.HaloDir, 'input_StarMaker.dat')
        ffile = os.path.join(self.HaloDir, 'inputfiles_StarMaker.dat')
        runcmd = OZYPATH + '/HaloFinder/program/GalFinder ./ > runGals.log'
        super().setup_HaloFinderRun(run=run, nvoisins=nvoisins, rhot=rhot, nhop=nhop, npart=npart,
                                    fudgepsilon=fudgepsilon, alphap=alphap, method=method, cenmethod='com',
                                    inputFile=inputFile, ffile=ffile, runCmd=runcmd, omega_b=True)