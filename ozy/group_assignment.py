import numpy as np
import numpy.ma as ma

def galaxies_to_halos(obj):
    """Assign galaxies to halos.
    
    This function finds the closest halo to a particular galaxy.
    """
    if 'linking_var' in obj._args:
        linking_variable = obj.linking_var
    else:
        linking_variable = 'total'
    
    if not obj._has_galaxies:
        # Check that the OZY object includes galaxies.
        return
    # Loop over all galaxies in catalogue.
    for galaxy in obj.galaxies:
        galaxy.parent_halo_index = -1
        distances                = np.full(obj.nhalos, np.infty)
        masses                   = np.zeros(obj.nhalos)
        # Compute configuration-space distance of the galaxy to all halos.
        for i in range(0, obj.nhalos):
            halo           = obj.halos[i]
            vec_dist       = galaxy.position - halo.position
            d              = np.linalg.norm(vec_dist)
            # TODO: Change to virial radius
            halo_linking_l = halo.radius[linking_variable]
            # halo_linking_l = 0.2 * halo.virial_quantities['radius']
            # If the galaxy is below the halo linking length, save distance.
            if d <= halo_linking_l:
                distances[i] = d
                masses[i]    = halo.virial_quantities['mass']
        # Assign most massive halo to galaxy.
        if not np.all(distances == np.infty):
            #galaxy.parent_halo_index = np.argmin(distances)
            galaxy.parent_halo_index = np.argmax(masses)
    
    for halo in obj.halos:
        halo.galaxy_index_list = []
    
    # Append galaxy indexes to each halo's galaxy list.
    for i in range(0, obj.ngalaxies):
        galaxy = obj.galaxies[i]
        if galaxy.parent_halo_index > -1:
            obj.halos[galaxy.parent_halo_index].galaxy_index_list.append(i)
    # Find lonely halos (i.e. those without a galaxy assigned).
    len_list = np.array([len(halo.galaxy_index_list) for halo in obj.halos])
    lonely_halos   = ma.masked_where(len_list > 0, obj.halos)
    if type(lonely_halos.mask) != np.ndarray:
        n_lonely_halos = len(obj.halos)
        lonely_halos = np.array(obj.halos)
        my_mask = np.full(n_lonely_halos,False)
    else:
        n_lonely_halos = len(lonely_halos[~lonely_halos.mask])
        my_mask = lonely_halos.mask
    if n_lonely_halos > 0:
        lonely_halos_index_list = [k for k in range(0, len(lonely_halos)) if lonely_halos[k] is not ma.masked]
        for i in range(0, len(obj.galaxies)):
            galaxy = obj.galaxies[i]
            if galaxy.parent_halo_index == -1:
                distances = np.full(n_lonely_halos, np.infty)
                for j in range(0, n_lonely_halos):
                    halo           = lonely_halos[~my_mask][j]
                    vec_dist       = galaxy.position - halo.position
                    d              = np.linalg.norm(vec_dist)
                    halo_linking_l = 1.0 * halo.radius[linking_variable]
                    # If the galaxy is below the halo linking length, save distance.
                    if d <= halo_linking_l:
                        distances[j] = d
                # Assign closest lonely halo to the galaxy
                if not np.all(distances == np.infty):
                    galaxy.parent_halo_index = lonely_halos_index_list[np.argmin(distances)]
                    obj.halos[galaxy.parent_halo_index].galaxy_index_list.append(i)


                
def central_galaxies(obj, central_mass_definition='total'):
    """Assign central galaxies.
    
    This function is in charge of determining the central galaxy as the most massive galaxy in a DM halo,
    and giving the status of sattellite galaxies to the rest.
    
    """
    if not obj._has_galaxies:
        return
    
    obj.central_galaxies   = []
    obj.satellite_galaxies = []
    
    for halo in obj.halos:
        if len(halo.galaxy_index_list) == 0:
            # This halo has no galaxies associated to it.
            continue
        galaxy_masses                                               = np.array([s.mass[central_mass_definition] for s in halo.galaxies])
        central_index                                               = np.argmax(galaxy_masses)
        obj.galaxies[halo.galaxy_index_list[central_index]].central = True