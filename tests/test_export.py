import numpy as np
import ozy
from ozy.export import unigrid_amr,basicexport2txt

obj = ozy.load('test_00035.hdf5')
virial_mass = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(virial_mass)
gal = obj.galaxies[progind]

r = 0.01 * obj.halos[gal.parent_halo_index].virial_quantities['radius'].to('code_length').d
pos = gal.position.to('code_length').d

grid = unigrid_amr(obj, group=gal, lmax=14, xmin=(pos[0]-r,'code_length'), xmax=(pos[0]+r,'code_length'),
                    ymin=(pos[1]-r,'code_length'), ymax=(pos[1]+r,'code_length'), zmin=(pos[2]-r,'code_length'),
                    zmax=(pos[2]+r,'code_length'))

print(grid.shape, grid.min(), grid.max())

basicexport2txt(obj,gal,rmax=(r,'code_length'),h=(0,'pc'),smoothmethod='level')

