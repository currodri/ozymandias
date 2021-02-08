import numpy as np
from pprint import pprint

grouptypes = dict(
    halo='halos',
    galaxy='galaxies',
    cloud='clouds'
)

class Group(object):
    """This is the parent class from which the rest of groups are derived."""
    
    def __init__(self, obj):
        """Initialise basic properties of group object."""
        self.obj      = obj
        self.ID       = -1
        self.npart    = 0
        self.units    = 'code'
        self.position = np.array([-1, -1, -1])
        self.radius   = -1
        self.shape = np.array([0, 0, 0])
        self.mass     = -1
        self.velocity = np.array([0, 0, 0])
        self.angular_mom = np.array([0, 0, 0])
        self.virial_mass = -1
        self.virial_radius = -1
        self.virial_temperature = -1
        self.virial_cvel = -1
        self.energies = np.array([0, 0, 0])
    def _code2physical(self, unit_l, unit_d, unit_t):
        """Convert halo properties from code to physical units."""
        if self.units == 'code':
            cm2kpc        = 3.24077929e-22 # Convert cm to kpc
            g2Msun        = 5e-34          # Convert g to solar masses
            unit_mass     = ((unit_d * unit_d) * unit_l) * unit_l
            # Update values.
            self.radius   = self.radius    * unit_l    *  cm2kpc
            self.mass     = self.mass      * unit_mass *  g2Msun
            self.position = self.position  * unit_l    *  cm2kpc
            self.velocity = self.velocity  * (unit_l/unit_t) *  1e-5 #  Express it in km/s
            self.virial_mass = self.virial_mass * unit_mass * g2Msun 
            self.virial_radius = self.virial_radius  * unit_l    *  cm2kpc
            # TODO: Is it virial temperature or circular velocity?
            #self.virial_temperature = 

            self.units    = 'physical'     # Update so we do not need to convert them again
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
    def __init__(self, obj):
        super(Cloud, self).__init__(obj)
        self.central = False
        self.galaxy  = None
        self.halo    = None

class Galaxy(Group):
    """Galaxy class which has the central boolean."""
    obj_type = 'galaxy'
    def __init__(self, obj):
        super(Galaxy, self).__init__(obj)
        self.central = False
        self.halo    = None
        self.clouds  = []
        self.cloud_index_list = []

class Halo(Group):
    """Halo class which has different levels of the halo hierarchy."""
    obj_type = 'halo'
    def __init__(self, obj):
        super(Halo, self).__init__(obj)
        self.spin = 0
        self.level              = -1
        self.host               = -1
        self.hostsub            = -1
        self.nsub               = -1
        self.nextsub            = -1
        self.galaxies           = []
        self.central_galaxy     = None
        self.satellite_galaxies = []
        self.galaxy_index_list = []

def create_new_group(obj, grouptype):
    """Simple function to create a new instance of a specified :class:`group.Group`.
    """
    if grouptype == 'halo':
        return Halo(obj)
    elif grouptype == 'galaxy':
        return Galaxy(obj)
    elif grouptype == 'cloud':
        return Cloud(obj)