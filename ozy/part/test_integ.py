import numpy as np
import sys
sys.path.append('/mnt/extraspace/currodri/Codes/ozymandias/ozy/part')
import part2
repository = '/mnt/extraspace/currodri/NUT/cosmoNUThd/output_00040'

varIDs = part2.io_ramses.hydroID()
part2.io_ramses.read_hydrofile_descriptor(repository,varIDs)

print(varIDs)

# reg = part2.geometrical_regions.region()
# reg.name = 'sphere'
# centre = part2.vectors.vector()
# centre.x,centre.y,centre.z = 0.70186129, 0.33655524, 0.34274765
# reg.centre = centre
# axis = part2.vectors.vector()
# axis.x,axis.y,axis.z = -0.14795277,  0.7482291 ,  0.6467327
# reg.axis = axis
# bulk = part2.vectors.vector()
# bulk.x,bulk.y,bulk.z = 0.00797963, -0.00414802, -0.00708832
# reg.bulk_velocity = bulk
# reg.rmin = 0.0
# reg.rmax = 0.0028789570555090905

# filt = part2.filtering.filter()

# glob_attrs = part2.part_integrator.part_region_attrs()
# glob_attrs.nvars = 4
# glob_attrs.nwvars = 1
# part2.part_integrator.allocate_part_regions_attrs(glob_attrs)
# glob_attrs.varnames.T.view('S128')[0] = b'star/mass'.ljust(128)
# glob_attrs.varnames.T.view('S128')[1] = b'star/sfr_10'.ljust(128)
# glob_attrs.varnames.T.view('S128')[2] = b'star/sfr_100'.ljust(128)
# glob_attrs.varnames.T.view('S128')[3] = b'dm/mass'.ljust(128)
# glob_attrs.wvarnames.T.view('S128')[0] = b'cumulative'.ljust(128)

# part2.part_integrator.integrate_region(repository,reg,filt,glob_attrs)

# print(glob_attrs.nstar,glob_attrs.ndm)
# print(glob_attrs.data)