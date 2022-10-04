import os
import numpy as np
import ozy
from ozy.projections import do_projection,plot_single_galaxy_projection,do_healpix_projection,plot_single_galaxy_healpix


obj = ozy.load('test_00035.hdf5')
virial_mass = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(virial_mass)
gal = obj.galaxies[progind]

# proj = do_projection(gal,['gas/density','gas/temperature','gas/v_cyl_z','gas/magnetic_magnitude','gas/streaming_heating','gas/grad_crp','gas/grav_crpfz'],window=(5,'kpc'),pov='z')
# proj.save_FITS('NUT_00035_z.fits')
# plot_single_galaxy_projection('NUT_00035_z.fits',['gas/density','gas/temperature','gas/v_cyl_z','gas/magnetic_magnitude','gas/streaming_heating','gas/grad_crp','gas/grav_crpfz'])

# proj = do_projection(gal,['gas/density','gas/temperature','gas/v_cyl_z','gas/v_sphere_r','gas/grav_crpfz','gas/grav_crpfrsphere','gas/grad_crpz','gas/grav_therpfz'],window=(0.1,'rvir'),pov='edgeon',weight=['gas/volume','star/cumulative'])
# proj.save_FITS('NUT_00035_edgeon.fits')
# plot_single_galaxy_projection('NUT_00035_edgeon.fits',['gas/density','gas/temperature','gas/v_cyl_z','gas/v_sphere_r','gas/grav_crpfz','gas/grav_crpfrsphere','gas/grad_crpz','gas/grav_therpfz'])

proj = do_projection(gal,['gas/density','gas/temperature','gas/eff_FKmag','gas/magnetic_magnitude','gas/streaming_heating','gas/grad_crp','gas/grav_crpfz'],window=(5,'kpc'),pov='faceon')
proj.save_FITS('NUT_00035_faceon.fits')
plot_single_galaxy_projection('NUT_00035_faceon.fits',['gas/density','gas/temperature','gas/eff_FKmag','gas/magnetic_magnitude','gas/streaming_heating','gas/grad_crp','gas/grav_crpfz'])

# proj = do_healpix_projection(gal,['gas/v_sphere_r','gas/density','gas/temperature','gas/metallicity'],r=(0.2,'rvir'),dr=(0.01,'rvir'))
# proj.save_FITS('NUT_00035_mollweide.fits')
# plot_single_galaxy_healpix('NUT_00035_mollweide.fits',['gas/v_sphere_r','gas/density','gas/temperature','gas/metallicity'])


# import yt
# import numpy as np
# import os, sys
# from subprocess import call
# import h5py
# from yt.fields.field_detector import FieldDetector

# def add_gradient_field(ds, gradient_field, direction):
#     '''Add a gradient field
#         Parameters
#         ----------
#         ds : Dataset
#         gradient_field : tuple
#         The field to compute the gradient of (ftype, fname), e.g. ("gas", "pressure")
#         direction : str
#         The gradient direction (x, y or z)
#         Returns
#         -------
#         fname : tuple
#         The name of the gradient field as a tuple.
#         '''
#     DIRECTION = 'xyz'.index(direction)
#     def _gradient_at_point(field, data):
#         # Get the position of the particles
#         pos = yt.ustack([data['index', k] for k in 'xyz'], axis=1)
#         Npart = pos.shape[0]
#         pos_neigh = yt.ustack((pos, pos), axis=0)
        
#         if isinstance(data, FieldDetector):
#             return np.zeros(Npart)
        
#         dx = data['index', 'dx'].to('code_length')
        
#         # Compute center of neighboring cells
#         pos_neigh[0, :, DIRECTION] -= dx
#         pos_neigh[1, :, DIRECTION] += dx
#         lvl = data['index', 'grid_level']
        
#         Npart = pos.shape[0]
#         tmp = np.zeros(Npart)
        
#         ILEFT = 0
#         IRIGHT = 1
        
#         ret = np.zeros((2, Npart))
#         remaining = np.ones((2, Npart), dtype=bool)
#         Nremaining = [Npart, Npart]
        
#         for i, subset in enumerate(data._current_chunk.objs):
#             mesh_data = subset[gradient_field].T.reshape(-1)
#             # Get cell index on left side & right side
#             for idir in (ILEFT, IRIGHT):
#                 if Nremaining[idir] == 0:
#                     continue
#                 tmp[:Nremaining[idir]] = subset.mesh_sampling_particle_field(
#                                 pos_neigh[idir, remaining[idir]].copy(),
#                                 mesh_data, lvlmax=lvl[remaining[idir]])
                    
#                 ret[idir, remaining[idir]] = tmp[:Nremaining[idir]]
                                                                             
#                 remaining[idir, remaining[idir]] = np.isnan(tmp[:Nremaining[idir]])
#                 Nremaining[idir] = remaining[idir].sum()
        
#         return data.ds.arr((ret[IRIGHT, :] - ret[ILEFT, :]) / dx / 2,
#                            mesh_data.units / dx.units)
#     ds.fields
#     units = '%s/%s' % (ds.field_info[gradient_field].units,
#                    ds.field_info['index', 'dx'].units)
    
#     ptype, fname = gradient_field
#     new_field = (ptype, f'{fname}_gradient_{direction}')
    
#     ds.add_field(new_field, function=_gradient_at_point,
#                 units=units, sampling_type='cell', force_override=True)
                 
#     return new_field

# def _Ptot(field,d):
#     Ptot = d['Pcr']+d['Pth']
#     return Ptot

# trans_matrix = np.array([[-0.33188517,  0.65337076, -0.68041082],
#                         [-0.89157098, -0.45288098,  0.        ],
#                         [-0.30814512,  0.60663455,  0.73283089]])
# trans_matrix = np.transpose(trans_matrix)
# # def _Z(field,d):
# #     Z = d['z'] - ds.arr(840.10831394, 'kpc')
# #     #Z = np.dot(Z,trans_matrix)
# #     return Z

# filename = '/mnt/extraspace/currodri/NUT/cosmoNUTcrmhd/output_00035/info_00035.txt'
# ds = yt.load(filename,fields=['Density','x-velocity','y-velocity','z-velocity','B_left_x','B_left_y','B_left,z','B_right_x','B_right_y','B_right_z','Pcr','Pth','Metallicity'])
# d=ds.all_data()


# ds.fields.ramses.Pcr.units='code_pressure'
# ds.fields.ramses.Pcr.output_units='code_pressure'
# ds.fields.ramses.Pth.units='code_pressure'
# ds.fields.ramses.Pth.output_units='code_pressure'

# #ds.add_field(("ramses","Z"),function=_Z,units="kpc")

# add_gradient_field(ds,('ramses','Pcr'),'z')
# add_gradient_field(ds,('ramses','Pth'),'z')

# Nbins=50
# height=5
# ZLOG=True
# PLOG=False
# disk=ds.disk(ds.domain_center+ds.arr(height,'kpc'),(0,0,1),(10,'kpc'),(height,'kpc'))
# dd=d.cut_region(["(np.isfinite(obj['Pcr_gradient_z']))"])
# regionPcr=ds.intersection((disk,dd))

# prj = yt.OffAxisProjectionPlot(ds,center=([1699.26980153,  823.10387361,  840.10831394], 'kpc'),
#                                     normal=[0.95133936,0.19649294,0.23736878],fields=[('ramses','Pcr_gradient_z')],width=(10,'kpc'), weight_field=("ramses", "cell_volume"))
# prj.save()