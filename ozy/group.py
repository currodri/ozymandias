import numpy as np
from pprint import pprint

grouptypes = dict(
    halo='halos',
    galaxy='galaxies',
    cloud='clouds'
)

info_blacklist = [
    'obj', 'halo', 'galaxies','clouds', 'satellites',
    'galaxy_index_list_end', 'galaxy_index_list_start','cloud_index_list_end','cloud_index_list_start']

# class GroupList(object):
#     """Class to hold particle index lists.
#     """
#     def __init__(self, name):
#         self.name = name

#     def __get__(self, instance, owner):
#         if not hasattr(instance, '_%s' % self.name) or isinstance(getattr(instance, '_%s' % self.name), int):
#             from ozy.loader import restore_single_list
#             restore_single_list(instance.obj, instance, self.name)
#         return getattr(instance, instance, self.name)
    
#     def __set__(self, instance, value):
#         setattr(instance, '_%s' % self.name, value)

class Group(object):
    """This is the parent class from which the rest of groups are derived."""

    # slist = GroupList('slist')
    
    def __init__(self, obj):
        """Initialise basic properties of group object."""
        self.obj      = obj
        self.ID       = -1
        self.npart    = 0
        self.units    = 'code'
        self.position = np.array([-1, -1, -1])
        self.radius   = {}
        self.shape = {}
        self.mass     = {}
        self.velocity = np.array([0, 0, 0])
        self.angular_mom = {}
        self.virial_quantities = {}
        self.energies = {}
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