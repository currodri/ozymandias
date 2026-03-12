import numpy as np
import ozy
from ozy.export import part2list

# Adjust this path if needed.
obj = ozy.load('test_00010.hdf5')
gal = obj.most_massive_system

# Use a small spherical region around the galaxy for a fast smoke test.
r = 0.01 * obj.halos[gal.parent_halo_index].virial_quantities['radius'].to('code_length').d

vars_to_extract = ['x', 'y', 'z', 'vx', 'vy', 'vz', 'mass', 'age', 'metallicity']

data = part2list(
    obj,
    group=gal,
    rmax=(r, 'code_length'),
    vars=vars_to_extract,
)

print('Requested vars:', vars_to_extract)
print('Output shape (nvars, nparticles):', data.shape)

if data.shape[1] > 0:
    print('First particle:')
    for i, vname in enumerate(vars_to_extract):
        print('  %s = %e' % (vname, data[i, 0]))
else:
    print('No particles selected. Try increasing rmax or changing the filter.')
