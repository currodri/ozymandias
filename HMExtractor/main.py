import numpy as np

from ozy.sim_attributes import SimulationAttributes

class OZY(object):
    """Master OZY class.
    OZY objects contain all the necessary references to halos, galaxies,
    and gas clouds in an individual simulation snapshot.

    It can be saved as a portable, standalone HDF5 file which allows
    general analysis without requiring the original snapshot.
    """
    def __init__(self, ds = 0, *args, *kwargs):
        self._args   = args
        self._kwargs = kwargs
        self._dtman  = 0
        self._ds     = 0

        self.units = dict(
            mass        = 'Msun',
            length      = 'kpccm',
            velocity    = 'km/s',
            time        = 'yr',
            temperature = 'K'
        )

        self.simulation  = SimulationAttributes()
        self.yt_dataset  = ds
        
        self.nhalos      = 0
        self.ngalaxies   = 0
        self.nclouds     = 0
        self.halos       = []
        self.galaxies    = []
        self.clouds      = []
        self.group_types = []

    @property
    def yt_dataset(self):
        """The yt dataset to perform actions on."""
        if isinstance(self._ds, int):
            raise Exception('No yt_dataset assigned!\nPlease assign '\
                            'one via `obj.yt_dataset=<YT DATASET>` ' \
                            'if you want to perform further analysis' \
                            'on the snapshot.')
        return self._ds
    @property
    def _has_galaxies(self):
        """Check if the dataset has galaxies."""
        if self.ngalaxies > 0:
            return True
        else:
            return False
    @property
    def _has_clouds(self):
        """Check if the dataset has gas clouds."""
        if self.nclouds > 0:
            return True
        else:
            return False
    
    @property
    def data_manager(self):
        """On demand DataManager class."""
        if isinstance(self.__dtman, int):
            from ozy.data_manager import DataManager
            self._dtman = DataManager(self)
        return self._dtman

    def _assign_simulation_attributes(self):
        """Assign simulation attributes to the OZY object, if it has not been done before."""
        self.simulation.assign_attributes(self)
    
    def _assign_groups(self):
        """Assign galaxies to halos and gas clouds to galaxies.
            Also connect halos with their central galaxy."""
        import ozy.group_assignment as assign
        assign.galaxies_to_halos(self)
        assign.central_galaxies(self)
        assign.clouds_to_galaxies(self)
    
    def _link_groups(self):
        """Two-way linking of objects."""
        from ozy.group_linking import link
        link.galaxies_to_halos(self)
        link.clouds_to_galaxies(self)
    
    def save(self, filename):
        """Save OZY object as HDF5 file."""
        from ozy.saver_tool import save
        save(self, filename)
    
    def build_HaloMaker(self, *args, *kwargs):
        """This is the central function of the OZY class for the HALOMAKER catalogues.

        This method is reponsible for:
        1) Calling the Fortran routines that cleans up the raw HaloMaker catalogues
        2) Creating halox, galaxies and gas clouds
        3) Linking objects through the chosen method
        4) Computing additional quantities
        5) Saving all as a clean HDF5 file

        """
        import ozy.group_assignment as assign
        import ozy.group_linking as link
        import ozy.read_HaloMaker as read_HM

        # TODO: Add the option to run HaloMaker if the brick files do not exist
        # import ozy.run_halomaker as run
        # run(self, 'halo')
        # run(self, 'galaxy')
        # run(self, 'cloud')

        # Read HaloMaker brick catalogues
        read_HM(self, 'halo')
        read_HM(self, 'galaxy')
        read_HM(self, 'cloud')

        # Make assignment
        assign.galaxies_to_halos(self)
        assign.clouds_to_galaxies(self)

        # Link objects between each other
        link.galaxies_to_halos(self)
        link.clouds_to_galaxies(self)

    def galaxies_summary(self, top=10):
        """Method to briefly print information for the most massive galaxies in the catalogue."""
        from ozy.tools import info_printer
        info_printer(self, 'galaxy', top)

    def halos_summary(self, top=10):
        """Method to briefly print information for the most massive halos in the catalogue."""
        from ozy.tools import info_printer
        info_printer(self, 'halo', top)

    def clouds_summary(self, top=10):
        """Method to briefly print information for the most massive clouds in the catalogue."""
        from ozy.tools import info_printer
        info_printer(self, 'cloud', top)