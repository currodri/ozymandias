import numpy as np
import ozy
from ozy.export import unigrid_amr

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

obj = ozy.load('test_00035.hdf5')
virial_mass = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(virial_mass)
gal = obj.galaxies[progind]

r = 0.1 * obj.halos[gal.parent_halo_index].virial_quantities['radius'].to('code_length').d
pos = gal.position.to('code_length').d

grid = unigrid_amr(obj, group=gal, lmax=14, xmin=(pos[0]-r,'code_length'), xmax=(pos[0]+r,'code_length'),
                    ymin=(pos[1]-r,'code_length'), ymax=(pos[1]+r,'code_length'), zmin=(pos[2]-r,'code_length'),
                    zmax=(pos[2]+r,'code_length'))

print(grid.shape, grid.min(), grid.max())

ax = plt.figure().add_subplot(projection='3d')
X, Y, Z = axes3d.get_test_data(0.05)

ax.contour(X, Y, Z, cmap=cm.coolwarm) 