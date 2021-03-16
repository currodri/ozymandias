import numpy as np
from pprint import pprint
from ozy.part2 import part2cube,part2sfr
from ozy.amr2 import amr2cube

MINIMUM_STARS_PER_GALAXY = 100
MINIMUM_DM_PER_HALO      = 0

grouptypes = dict(
    halo='halos',
    galaxy='galaxies',
)

info_blacklist = [
    'obj', 'halo', 'galaxies', 'satellites',
    'galaxy_index_list_end', 'galaxy_index_list_start']

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
    
    @property
    def _valid(self):
        """Check against the minimum number of particles to see if
        this object is 'valid'."""
        if self.obj_type == 'halo' and self.ndm < MINIMUM_DM_PER_HALO:
            return False
        elif self.obj_type == 'galaxy' and self.nstar < MINIMUM_STARS_PER_GALAXY:
            return False
        else:
            return True

    def _info(self):
        """Method to quickly print out object attributes."""
        attrdict = {}
        for k,v in six.iteritems(self.__dict__):
            if k in info_blacklist:
                continue
            attrdict[k] = v
        pprint(attrdict)
        attrdict = None

class Galaxy(Group):
    """Galaxy class which has the central boolean."""
    obj_type = 'galaxy'
    def __init__(self, obj):
        super(Galaxy, self).__init__(obj)
        self.central = False
        self.halo    = None
        self.sfr     = {}
    def _process_galaxy(self):
        """Process each galaxy after creation. This means
        calculating the total mass, and then calculate the rest of masses,
        radial quantities, velocity dispersions, angular momentum...
        """
        self._calculate_masses()
        self._calculate_velocity_dispersions()
        self._calculate_gas_quantities()
        self._calculate_star_quantities()
    def _calculate_masses(self):
        """Calculate various total masses"""
        output_path = self.obj.simulation.fullpath

        # Region for selection -- Using: Radius of gas+stars structure
        r_region = self.radius['total'].in_units('code_length').d
        x_center = self.position.in_units('code_length')[0].d
        y_center = self.position.in_units('code_length')[1].d
        z_center = self.position.in_units('code_length')[2].d

        # subroutine integratesphere(repository,ageweight,periodic,star,xcenter,ycenter,zcenter,radius)
        mass_dm, ndm = part2cube.integratesphere(output_path, False, True, False, x_center, 
                                                    y_center, z_center, r_region)
        mass_star, nstar = part2cube.integratesphere(output_path, False, True, True, x_center, 
                                                        y_center, z_center, r_region)
        # subroutine integratesphere(repository,namevar,xcenter,ycenter,zcenter,radius,lmax)
        mass_gas = amr2cube.integratesphere(output_path, 'mass', x_center, y_center, z_center,
                                            r_region, 10)
        mass_baryon = mass_gas + mass_star

        self.mass['dm'] = self.obj.yt_dataset.quan(mass_dm, 'code_mass')
        self.mass['gas']     = self.obj.yt_dataset.quan(mass_gas, 'code_mass')
        self.mass['stellar'] = self.obj.yt_dataset.quan(mass_star, 'code_mass')
        self.mass['baryon']  = self.obj.yt_dataset.quan(mass_baryon, 'code_mass')

        self.ndm = ndm
        self.nstar = nstar

        self.gas_fraction = 0.0
        if self.mass['baryon'] > 0:
            self.gas_fraction = self.mass['gas'].d / self.mass['baryon'].d
            
    def _calculate_gas_quantities(self):
        """Compute gas quantities: Metallicity, Temperature..."""
        # TODO
        return
    def _calculate_star_quantities(self):
        """Calculate star quantities..."""
        # TODO
        indicators = np.array([0.01, 0.1]) # in Gyr. These should be taken as arguments in the future
        output_path = self.obj.simulation.fullpath

        # Region for selection -- Using: Radius of gas+stars structure
        r_region = self.radius['total'].in_units('code_length').d
        x_center = self.position.in_units('code_length')[0].d
        y_center = self.position.in_units('code_length')[1].d
        z_center = self.position.in_units('code_length')[2].d
        
        sfr = part2sfr.sphere(output_path, x_center, y_center,
                                z_center, r_region, indicators, 
                                len(indicators))

        for i in range(0, len(indicators)):
            self.sfr[str(int(indicators[i]*1e+3))+'Myr'] = self.obj.yt_dataset.quan(sfr[i],'Msun/yr')

    def _calculate_velocity_dispersions(self):
        """Calculate velocity dispersions for the various components."""
        # TODO
        return
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