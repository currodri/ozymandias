import numpy as np
import amr2prof

l = amr2prof.vectors.vector()
mat = np.zeros((3,3),order='F')
errormsg = 0
l.x,l.y,l.z = 0.5773502691896258,0.5773502691896258,0.5773502691896258
amr2prof.coordinate_systems.new_z_coordinates(l,mat,errormsg)
b = np.array([0.5773502691896258,0.5773502691896258,0.5773502691896258], order='F')
# b = np.array([1/np.sqrt(2),1/np.sqrt(2),0],order='F')
print('np.pi/4= ',np.pi/4)
print(mat)
print(mat.dot(b))