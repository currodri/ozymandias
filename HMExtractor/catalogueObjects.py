import numpy as np
import six
from pprint import pprint

category_mapper = dict(
    mass = 'masses',
    radius = 'radii',
    sigma = 'velocity_dispersion',
    metallicity = 'temperatures'
)

class GroupProperty(object):
    """Class to return default values for the quantities held in the category_mapper dictionary."""
    def __init__(self, source_dict, name):
        self.source_dict = source_dict
        self.name = name
    def __get__(self, instance, owner):
        pass

class Group(object):
    """This is the parent class from which the rest of groups are derived."""
    def __init__(self):
        """Initialise basic properties of group object."""
        self.ID = -1
        self.units = 'code'
        self.position = np.array([-1, -1, -1])
        self.radius = -1
        self.mass = -1
        self.velocity = np.array([-1, -1, -1])
        self.angular_mom = np.array([-1, -1, -1])
    def _code2physical(self, unit_l, unit_d, unit_t):
        """Convert halo properties from code to physical units."""
        if self.units == 'code':
            cm2kpc=3.24077929e-22 # Convert cm to kpc
            g2Msun=5e-34 # Convert g to solar masses
            unit_mass = ((unit_d * unit_d) * unit_l) * unit_l
            # Update values
            self.radius = self.radius * unit_l * cm2kpc
            self.mass = self.mass * unit_mass * g2Msun
            self.position = self.position * unit_l * cm2kpc
            self.velocity = self.velocity * (unit_l/unit_t) * 1e-5 # Express it in km/s

            self.units = 'physical' # Update so we do not need to convert them again
    def _info(self):
        """Method to quickly print out object attributes."""
        attrdict = {}
        for k,v in six.iteritems(self.__dict__):
            if k in info_blacklist:
                continue
            attrdict[k] = v
        pprint(attrdict)
        attrdict = None

    
class Cloud(Group):
    """Group made of gas particles from R2G."""
    obj_type = 'cloud'

class Galaxy(Group):
    """Galaxy class which has the central boolean."""
    obj_type = 'galaxy'
    def __init__(self):
        self.central = False
        self.halo = None

class Halo(Group):
    """Halo class which has different levels of the halo hierarchy."""
    obj_type = 'halo'
    def __init__(self):
        self.galaxies = []
        self.central_galaxy = None
        self.satellite_galaxies = []

