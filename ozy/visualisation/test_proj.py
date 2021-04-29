import numpy as np
import sys
sys.path.append('/mnt/extraspace/currodri/Codes/ozymandias/ozy/visualisation')
import projections

repository = '/mnt/extraspace/currodri/NUT/cosmoNUThd/output_00044'
centre = projections.vectors.vector()
centre.x,centre.y,centre.z = 0.70336995, 0.33568012, 0.341516
los_axis = projections.vectors.vector()
los_axis.x,los_axis.y,los_axis.z = -0.35676077,  0.40252785,  0.84302614
up_vector = projections.vectors.vector()
region_size = np.array([2*0.0028299868106842,2*0.0028299868106842],order='F')
distance = 0.0028299868106842
far_cut_depth = 0.0028299868106842
map_max_size = 1024

bulk = projections.vectors.vector()
bulk.x,bulk.y,bulk.z = 28.26907349, -20.34032822, -25.46469116
cam = projections.obs_instruments.init_camera(centre,los_axis,up_vector,region_size,distance,far_cut_depth,map_max_size)

trans = np.zeros((3,3),order='F')
projections.obs_instruments.los_transformation(cam,trans)
axis = np.array([-0.35676077,  0.40252785,  0.84302614])
print(trans, np.dot(trans,axis))
print(np.dot(trans.T,[0,0,1]))

# proj = projections.amr_map.projection_handler()
# proj.type = 'faceon'
# proj.nvars = 3
# proj.weightvar = 'density'
# projections.amr_map.allocate_projection_handler(proj)
# proj.varnames.T.view('S128')[0] = b'density'.ljust(128)
# proj.varnames.T.view('S128')[1] = b'temperature'.ljust(128)
# proj.varnames.T.view('S128')[2] = b'metallicity'.ljust(128)


# projections.amr_map.projection(repository,cam,bulk,proj)

