import numpy as np
import amr2prof

repository = '/mnt/extraspace/currodri/NUT/cosmoNUThd/output_00044'
reg = amr2prof.geometrical_regions.region()
reg.name = 'cylinder'
centre = amr2prof.vectors.vector()
centre.x,centre.y,centre.z = 0.70336995, 0.33568012, 0.341516
reg.centre = centre
axis = amr2prof.vectors.vector()
axis.x,axis.y,axis.z = -0.35676077,  0.40252785,  0.84302614
reg.axis = axis
bulk = amr2prof.vectors.vector()
bulk.x,bulk.y,bulk.z = 28.26907349, -20.34032822, -25.46469116
reg.bulk_velocity = bulk
reg.rmin = 0.0
reg.rmax = 0.0028299868106842
reg.zmin = -0.0028299868106842
reg.zmax = 0.0028299868106842
lim = np.zeros((3,2),order='F')
amr2prof.geometrical_regions.limits(reg,lim)
print(lim)

box_cyl = np.zeros((8,3))
axis = np.array([-0.35676077,  0.40252785,  0.84302614])
centre = np.array([0.70336995, 0.33568012, 0.341516])
v = np.array([0.156876702141575,0.915407738832915,-0.370699840855096])
w = np.array([0.920929667244383,0,0.389728813393556])
for i in range(0,3):
    box_cyl[0,i] = centre[i] + 0.0028299868106842*axis[i] + np.sqrt(2)*0.0028299868106842*v[i]
    box_cyl[1,i] = centre[i] + 0.0028299868106842*axis[i] - np.sqrt(2)*0.0028299868106842*v[i]
    box_cyl[2,i] = centre[i] + 0.0028299868106842*axis[i] + np.sqrt(2)*0.0028299868106842*w[i]
    box_cyl[3,i] = centre[i] + 0.0028299868106842*axis[i] - np.sqrt(2)*0.0028299868106842*w[i]
    box_cyl[4,i] = centre[i] - 0.0028299868106842*axis[i] + np.sqrt(2)*0.0028299868106842*v[i]
    box_cyl[5,i] = centre[i] - 0.0028299868106842*axis[i] - np.sqrt(2)*0.0028299868106842*v[i]
    box_cyl[6,i] = centre[i] - 0.0028299868106842*axis[i] + np.sqrt(2)*0.0028299868106842*w[i]
    box_cyl[7,i] = centre[i] - 0.0028299868106842*axis[i] - np.sqrt(2)*0.0028299868106842*w[i]
print(box_cyl)
print(np.max(box_cyl[:,0]),np.max(box_cyl[:,1]),np.max(box_cyl[:,2]))
print(np.min(box_cyl[:,0]),np.min(box_cyl[:,1]),np.min(box_cyl[:,2]))
amr = amr2prof.io_ramses.amr_info()
sim = amr2prof.io_ramses.sim_info()

amr2prof.io_ramses.init_amr_read(repository,amr,sim)
amr.lmax =19
amr2prof.io_ramses.get_cpu_map(reg,amr)

trans_matrix = np.zeros((3,3),order='F')
errormsg =0
amr2prof.coordinate_systems.new_z_coordinates(reg.axis,trans_matrix,errormsg)
print(trans_matrix)
print(np.dot(trans_matrix,np.array([-0.35676077,  0.40252785,  0.84302614], order='F')))

