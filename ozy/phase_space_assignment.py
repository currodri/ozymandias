import numpy as np
import numpy.ma as ma

def galaxies_to_halos(obj):
    """
    Assign galaxies to halos following a phase-space approach.

    This functions finds the closest DM halo to a galaxy group in
    phase space.
    """
    if not obj._has_galaxies:
        # Check that the OZY object includes galaxies.
        return
    
    # Loop over all galaxies in catalogue
    for i in range(0, obj.ngalaxies):
        galaxy = obj.galaxies[i]
        gal_linked = False
        # 1) Get an estimate of the interparticle distance
        galaxy.parent_halo_index = -1
        V                        = (4./3.)*np.pi*galaxy.radius['200crit'].to('code_length')
        n                        = galaxy.npart / V
        rmean                    = 0.2 / (n**(1./3.))
        sigV                     = galaxy.velocity['sigV'].to('km/s')
        x_gal                    = galaxy.position['COM'].to('code_length')
        v_gal                    = galaxy.velocity['COM'].to('km/s')
        distances                = np.full(obj.nhalos, np.infty)
        # 2) Compute phase-space metric for all halos
        for j in range(0, obj.nhalos):
            halo   = obj.halos[j]
            x_dist = x_gal - halo.position['COM'].to('code_length')
            x_dist = np.linalg.norm(x_dist)
            v_dist = v_gal - halo.velocity['COM'].to('km/s')
            v_dist = np.linalg.norm(v_dist)
            D2     = (x_dist**2.) / (rmean**2.) + (v_dist**2.) / (sigV**2.)
            distances[j] = D2
        # 3) Find the closest one that is within 1 rvir of the DM halo
        while not gal_linked:
            if np.all(distances == np.infty):
                galaxy.parent_halo_index = ind
                gal_linked = True
                break
            ind = np.argmin(distances)
            halo   = obj.halos[ind]
            x_dist = x_gal - halo.position['COM'].to('code_length')
            d = np.linalg.norm(x_dist)
            if d <= halo.radius['200crit'].to('code_length').d:
                galaxy.parent_halo_index = ind
                gal_linked = True
            else:
                distances[ind] = np.infty
        # 4) Append galaxy indexes to each halo's galaxy list
        obj.halos[galaxy.parent_halo_index].galaxy_index_list.append(i)

def central_galaxies(obj, central_mass_definition='200crit'):
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