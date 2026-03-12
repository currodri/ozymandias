import numpy as np
import numpy.ma as ma
from scipy.spatial import cKDTree


def _reset_halo_assignments(obj):
    for halo in obj.halos:
        halo.central_galaxy = None
        halo.satellite_galaxies = []
        halo.galaxy_index_list = []

    for galaxy in obj.galaxies:
        galaxy.halo = None
        galaxy.parent_halo_index = -1


def _validate_halo_assignments(obj):
    central_counts = np.zeros(obj.nhalos, dtype=int)

    for galaxy in obj.galaxies:
        if galaxy.parent_halo_index < 0 or not galaxy.central:
            continue
        central_counts[galaxy.parent_halo_index] += 1

    duplicate_halos = np.where(central_counts > 1)[0]
    if duplicate_halos.size > 0:
        halo_ids = [obj.halos[i].ID for i in duplicate_halos[:10]]
        raise RuntimeError(
            'Multiple central galaxies were assigned to the same halo: '
            f'{halo_ids}'
        )

def satellites_to_galaxies(obj):
    """Assign satellite galaxies to their central galaxy.
    This function uses the hierarchy of GalFinder to link satellites
    to their central galaxy.
    """
    if not obj._has_galaxies:
        return

    for galaxy in obj.galaxies:
        satellites = galaxy.substructure_list
        galaxy.nsubs = len(satellites)
        galaxy.subs_list = satellites

def galaxies_to_halos(obj):

    halos     = obj.halos
    nhalos    = obj.nhalos
    galaxies  = obj.galaxies
    ngalaxies = obj.ngalaxies

    _reset_halo_assignments(obj)

    # 1. Setup mass arrays
    assigned = np.full(ngalaxies, False, dtype=bool)
    if obj.add_satellites:
        galMasses = np.array([gal._get_totmass_withsubs() for gal in galaxies])
    else:
        galMasses = np.array([gal.mass['total'] for gal in galaxies])
    haloMassesvir = np.array([halo.virial_quantities['mass'] for halo in halos])
    imassG = np.argsort(galMasses)[::-1]
    imassH = np.argsort(haloMassesvir)[::-1]

    # 2. Extract positions and radii as plain float arrays in Mpc (once, outside the loop)
    #    This avoids repeated unyt unit-conversion overhead inside the distance calls.
    boxsize  = float(obj.simulation.boxsize.to('Mpc'))
    gal_pos  = np.array([gal.position.to('Mpc').d  for gal in galaxies])   # (ngalaxies, 3)
    halo_pos = np.array([halo.position.to('Mpc').d for halo in halos])     # (nhalos, 3)
    rvirInside = 0.3 * np.array([halo.virial_quantities['radius'].to('Mpc').d for halo in halos])

    central_indices = np.array([i for i, gal in enumerate(galaxies) if gal.central], dtype=int)
    if central_indices.size == 0:
        central_indices = np.arange(ngalaxies, dtype=int)
        if obj._kwargs.get('verbose', False):
            print('No central galaxies flagged by HaloMaker; falling back to all galaxies for halo assignment.', flush=True)

    # 3. Mass rank: mass_rank[gi] = rank of galaxy gi in descending mass order (0 = most massive)
    mass_rank = np.empty(ngalaxies, dtype=int)
    mass_rank[imassG] = np.arange(ngalaxies)

    # 4. Build a KD-tree over galaxy positions with periodic boundary conditions.
    #    query_ball_point replaces the O(N_galaxies) inner Python loop, reducing
    #    total complexity from O(N_halos * N_galaxies) to O(N_galaxies*log + N_halos*log).
    tree = cKDTree(gal_pos[central_indices], boxsize=boxsize)

    # 5. Loop over halos from most to least massive
    for i in range(nhalos):
        hi   = imassH[i]
        halo = halos[hi]

        # Find all galaxies within 0.3 Rvir in one tree query (periodic BC handled by cKDTree)
        local_candidates = tree.query_ball_point(halo_pos[hi], rvirInside[hi])
        if not local_candidates:
            continue

        candidates = central_indices[local_candidates]

        # Filter to unassigned candidates, then pick the most massive
        unassigned = [gi for gi in candidates if not assigned[gi]]
        if not unassigned:
            continue

        indGal = min(unassigned, key=lambda gi: mass_rank[gi])
        assigned[indGal] = True
        galaxies[indGal].parent_halo_index = hi
        galaxies[indGal].halo = halo
        halo.central_galaxy = galaxies[indGal]
        halo.galaxy_index_list.append(indGal)

        # Propagate the central galaxy's satellites to the same halo and keep them
        # out of later central-galaxy assignment.
        if galaxies[indGal].nsubs > 0:
            for satellite in galaxies[indGal].subs_list:
                assigned[satellite._index] = True
                satellite.parent_halo_index = hi
                satellite.halo = halo
                halo.galaxy_index_list.append(satellite._index)
                halo.satellite_galaxies.append(satellite)

    _validate_halo_assignments(obj)

    # 6. Diagnostics
    if obj._kwargs.get('verbose', False):
        assigned_count      = np.sum(assigned)
        halos_with_galaxies = sum(1 for halo in halos if halo.central_galaxy is not None)
        print(f"Assigned {assigned_count}/{ngalaxies} galaxies to halos", flush=True)
        print(f"Assigned galaxies to {halos_with_galaxies}/{nhalos} halos", flush=True)

                
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
    obj.unassigned_galaxies = []
    obj.starless_halos = []
    
    # 1. Distribute galaxies in central, satellite and unassigned lists.
    for galaxy in obj.galaxies:
        if galaxy.halo is not None and galaxy.halo.central_galaxy is galaxy:
            obj.central_galaxies.append(galaxy)
        elif galaxy.halo is not None:
            obj.satellite_galaxies.append(galaxy)
        else:
            if not hasattr(obj, 'unassigned_galaxies'):
                obj.unassigned_galaxies = []
            obj.unassigned_galaxies.append(galaxy)

    # 2. Now get the halos which have no galaxies assigned
    for halo in obj.halos:
        if halo.central_galaxy is None:
            obj.starless_halos.append(halo)
