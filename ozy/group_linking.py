import numpy as np

def galaxies_to_halos(obj):
    """Link galaxies and halos to one another.
    
    This function creates the two-way link of galaxies and halos.
    It is run during the creating and loading of the OZY catalogue file.
    
    """
    
    if not obj._has_galaxies:
        return
    
    # First let's do halos.
    for halo in obj.halos:
        halo.galaxies = []
        for gal_index in halo.galaxy_index_list:
            halo.galaxies.append(obj.galaxies[gal_index])
    
    # Now do the same but for galaxies.
    for galaxy in obj.galaxies:
        if galaxy.parent_halo_index > -1:
            galaxy.halo = obj.halos[galaxy.parent_halo_index]
        else:
            continue

def clouds_to_galaxies(obj):
    """Link gas clouds and galaxies to one another.
    
    This function creates the two-way linking of gas clouds and galaxies.
    It is run during the creating and loading of the OZY catalogue file.
    """
    
    if not obj._has_clouds:
        return
    
    # First let's do galaxies.
    for galaxy in obj.galaxies:
        galaxy.clouds = []
        for cloud_index in galaxy.cloud_index_list:
            galaxy.clouds.append(obj.clouds[cloud_index])
    
    # Now the same but for gas clouds.
    for cloud in obj.clouds:
        if cloud.parent_galaxy_index > -1:
            cloud.galaxy = obj.galaxies[cloud.parent_galaxy_index]
        else:
            cloud.galaxy = None

def create_sublists(obj):
    """Create sublists of objects.
    
    This will create the following sublists in the OZY catalogue:
        - central_galaxies
        - satellite_galaxies
        - unassigned_galaxies (those without a halo)
    
    """
    
    if not obj._has_galaxies:
        return
    
    obj.central_galaxies = []
    obj.satellite_galaxies = []
    
    # Assign to each halo the ID of its central galaxy.
    for halo in obj.halos:
        halo.central_galaxy = -1
        for galaxy in halo.galaxies:
            if galaxy.central:
                halo.central_galaxy = galaxy.ID
    
    # Distribute galaxies in central, satellite and unassigned lists.
    for galaxy in obj.galaxies:
        if galaxy.central and galaxy.halo is not None:
            obj.central_galaxies.append(galaxy)
        elif galaxy.halo is not None:
            obj.satellite_galaxies.append(galaxy)
        else:
            if not hasattr(obj, 'unassigned_galaxies'):
                obj.unassigned_galaxies = []
            obj.unassigned_galaxies.append(galaxy)