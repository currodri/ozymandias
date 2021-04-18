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
print(cam.centre)
# bbox = projections.geometrical_regions.region()

# projections.obs_instruments.get_bounding_box(cam,bbox)
# print(bbox)
print(projections.obs_instruments.get_required_resolution(cam))

projections.amr_map.projection(repository,cam,bulk)

# repository = '/mnt/extraspace/currodri/NUT/cosmoNUThd/output_00044'
# reg = projections.geometrical_regions.region()
# reg.name = 'cube'
# centre = projections.vectors.vector()
# centre.x,centre.y,centre.z = 0.70336995, 0.33568012, 0.341516
# reg.centre = centre
# axis = projections.vectors.vector()
# axis.x,axis.y,axis.z = -0.35676077,  0.40252785,  0.84302614
# reg.axis = axis
# bulk = projections.vectors.vector()
# bulk.x,bulk.y,bulk.z = 28.26907349, -20.34032822, -25.46469116
# reg.bulk_velocity = bulk
# reg.xmin = 0.70336995 - 0.0028299868106842
# reg.xmax = 0.70336995 + 0.0028299868106842
# reg.ymin = 0.33568012 - 0.0028299868106842
# reg.ymax = 0.33568012 + 0.0028299868106842
# reg.zmin = 0.341516 - 0.0028299868106842
# reg.zmax = 0.341516 + 0.0028299868106842

# data = projections.amr_map.map_data()

# projections.amr_map.projection(repository,reg,data,0)