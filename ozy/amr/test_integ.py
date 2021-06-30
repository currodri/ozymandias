import numpy as np
import sys
sys.path.append('/mnt/extraspace/currodri/Codes/ozymandias/ozy/amr')
import amr2

repository = '/mnt/extraspace/currodri/NUT/cosmoNUThd/output_00040'

reg = amr2.geometrical_regions.region()
reg.name = 'sphere'
centre = amr2.vectors.vector()
centre.x,centre.y,centre.z = 0.70186129, 0.33655524, 0.34274765
reg.centre = centre
axis = amr2.vectors.vector()
axis.x,axis.y,axis.z = -0.14795277,  0.7482291 ,  0.6467327
reg.axis = axis
bulk = amr2.vectors.vector()
bulk.x,bulk.y,bulk.z = 0.00797963, -0.00414802, -0.00708832
reg.bulk_velocity = bulk
reg.rmin = 0.0
reg.rmax = 0.0028789570555090905

filt = amr2.filtering.filter()

glob_attrs = amr2.amr_integrator.region_attrs()
glob_attrs.nvars = 3
glob_attrs.nwvars = 1
amr2.amr_integrator.allocate_regions_attrs(glob_attrs)
glob_attrs.varnames.T.view('S128')[0] = b'mass'.ljust(128)
glob_attrs.varnames.T.view('S128')[1] = b'thermal_energy'.ljust(128)
glob_attrs.varnames.T.view('S128')[2] = b'ang_momentum'.ljust(128)
glob_attrs.wvarnames.T.view('S128')[0] = b'cumulative'.ljust(128)

amr2.amr_integrator.integrate_region(repository,reg,filt,glob_attrs)

print(glob_attrs.data)
