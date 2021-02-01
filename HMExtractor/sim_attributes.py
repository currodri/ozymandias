import numpy as np

class SimulationAttributes(object):
    """Class that contains the attributes of the simulation."""

    def __init__(self):
        pass
    
    def assign_attributes(self, obj):
        """Get attributes from yT and assign them to OZY object."""
        # TODO: Make sure to add an option that does not use yT but just
        # reads in the info_xxxxx.txt file.
        ds = obj.yt_dataset

        self.redshift        = ds.current_redshift
        self.scale_factor    = 1. / (1. + self.redshift)
        self.time            = ds.current_time
        self.omega_matter    = ds.omega_matter
        self.omega_lambda    = ds.omega_lambda
        self.fullpath        = ds.fullpath
        self.hubble_constant = ds.hubble_constant
        self.parameters      = ds.parameters
        self.boxsize = ds.domain_width[0].to(obj.units['length'])
        