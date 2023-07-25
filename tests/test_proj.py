import os
import numpy as np
import ozy
from ozy.projections import do_projection,plot_single_galaxy_projection,do_healpix_projection,plot_single_galaxy_healpix


obj = ozy.load('test_00035.hdf5')
virial_mass = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(virial_mass)
gal = obj.galaxies[progind]

# proj = do_projection(gal,['gas/density','gas/temperature','star/sdensity',
#                           'gas/neighbour_accuracy','gas/grad_crp','gas/eff_FKmag'],
#                             window=(5,'kpc'),pov='faceon',
#                             filter_conds=['none','star_age/<=/0.01/Myr'],
#                             filter_name=['all','young'],remove_subs=True)
# proj.save_FITS('NUT_00035_faceon.fits')
# plot_single_galaxy_projection('NUT_00035_faceon.fits',['gas/density','gas/temperature','star/sdensity',
#                                                   'gas/neighbour_accuracy','gas/grad_crp','gas/eff_FKmag'],
#                               filter_name='young')
# plot_single_galaxy_projection('NUT_00035_faceon.fits',['gas/density','gas/temperature','star/sdensity',
#                                                   'gas/neighbour_accuracy','gas/grad_crp','gas/eff_FKmag'],
#                               filter_name='all')

# proj = do_projection(gal,['gas/density','gas/temperature','star/sdensity','dm/sdensity'],
#                             window=(5,'kpc'),pov='edgeon')
# proj.save_FITS('NUT_00035_edgeon.fits')
# plot_single_galaxy_projection('NUT_00035_edgeon.fits',['gas/density','gas/temperature',
#                                                        'star/sdensity','dm/sdensity'])

# proj = do_projection(gal,['gas/density','gas/temperature','star/sdensity',
#                           'gas/neighbour_accuracy','gas/grad_crp','gas/eff_FKmag'],
#                             window=(5,'kpc'),pov='z',
#                             filter_conds=['none','star_age/<=/0.01/Myr'],
#                             filter_name=['all','young'])
# proj.save_FITS('NUT_00035_z.fits')
# plot_single_galaxy_projection('NUT_00035_z.fits',['gas/density','gas/temperature','star/sdensity',
#                                                   'gas/neighbour_accuracy','gas/grad_crp','gas/eff_FKmag'],
#                               filter_name='young')
# plot_single_galaxy_projection('NUT_00035_z.fits',['gas/density','gas/temperature','star/sdensity',
#                                                   'gas/neighbour_accuracy','gas/grad_crp','gas/eff_FKmag'],
#                               filter_name='all')

# proj = do_projection(gal,['gas/density','gas/temperature','star/sdensity',
#                           'gas/neighbour_accuracy','gas/grad_crp','gas/eff_FKmag'],
#                             window=(5,'kpc'),pov='x',
#                             filter_conds=['none','star_age/<=/0.01/Myr'],
#                             filter_name=['all','young'])
# proj.save_FITS('NUT_00035_x.fits')
# plot_single_galaxy_projection('NUT_00035_x.fits',['gas/density','gas/temperature','star/sdensity',
#                                                   'gas/neighbour_accuracy','gas/grad_crp','gas/eff_FKmag'],
#                               filter_name='young')
# plot_single_galaxy_projection('NUT_00035_x.fits',['gas/density','gas/temperature','star/sdensity',
#                                                   'gas/neighbour_accuracy','gas/grad_crp','gas/eff_FKmag'],
#                               filter_name='all')
# proj = do_projection(gal,['gas/density','gas/temperature','star/sdensity',
#                           'gas/magnetic_magnitude','gas/grad_crp','gas/eff_FKmag'],
#                             window=(5,'kpc'),pov='y')
# proj.save_FITS('NUT_00035_y.fits')
# plot_single_galaxy_projection('NUT_00035_y.fits',['gas/density','gas/temperature','star/sdensity',
#                                                   'gas/magnetic_magnitude','gas/grad_crp','gas/eff_FKmag'])
# plot_single_galaxy_projection('NUT_00035_y.fits',['gas/density','gas/temperature','star/sdensity',
#                                                   'gas/neighbour_accuracy','gas/grad_crp','gas/eff_FKmag'])
# proj = do_projection(gal,['gas/density','gas/temperature','gas/alfven_speed',
#                           'gas/diffusion_speed','gas/alfvendiff_ratio',
#                           'gas/grav_crpfrsphere','gas/sigma','gas/grav_therpfz'],
#                      window=(0.1,'rvir'),pov='edgeon',weight=['gas/density','star/cumulative'],
#                      type_projection='gauss_deposition',nexp_factor=0.4)
# proj.save_FITS('NUT_00035_edgeon.fits')
# plot_single_galaxy_projection('NUT_00035_edgeon.fits',['gas/density','gas/temperature','gas/alfven_speed',
#                                                        'gas/diffusion_speed','gas/alfvendiff_ratio',
#                                                        'gas/grav_crpfrsphere','gas/sigma',
#                                                        'gas/grav_therpfz'],type_scale='galaxy',smooth=False)

proj = do_projection(gal,['gas/density','gas/temperature','gas/alfven_speed',
                          'gas/cr_pressure','star/sfr_surface_100','dm/sdensity'],
                     window=(0.1,'rvir'),pov='edgeon',weight=['gas/density','star/cumulative'],
                     type_projection='gauss_deposition',nexp_factor=1.0,lmin=16,lmax=16)
proj.save_FITS('NUT_00035_edgeon.fits')
plot_single_galaxy_projection('NUT_00035_edgeon.fits',['gas/density','gas/temperature','gas/alfven_speed',
                          'gas/cr_pressure','star/sfr_surface_100','dm/sdensity'],type_scale='galaxy',smooth=False)

# proj = do_projection(gal,['gas/density','gas/temperature','gas/alfven_speed',
#                           'gas/metallicity','star/sfr_surface_100','dm/sdensity'],
#                      window=(0.1,'rvir'),pov='faceon',weight=['gas/density','star/cumulative'],
#                      type_projection='gauss_deposition',nexp_factor=2)
# proj.save_FITS('NUT_00035_faceon.fits')
# plot_single_galaxy_projection('NUT_00035_faceon.fits',['gas/density','gas/temperature','gas/alfven_speed',
#                           'gas/metallicity','star/sfr_surface_100','dm/sdensity'],type_scale='galaxy',smooth=False)

# proj = do_projection(gal,['gas/density','gas/temperature','gas/eff_FKmag','gas/magnetic_magnitude','gas/streaming_heating','gas/grad_crp'],window=(1,'kpc'),pov='faceon')
# proj.save_FITS('NUT_00035_faceon.fits')
# plot_single_galaxy_projection('NUT_00035_faceon.fits',['gas/density','gas/temperature','gas/eff_FKmag','gas/magnetic_magnitude','gas/streaming_heating','gas/grad_crp'])

# proj = do_healpix_projection(gal,['gas/massflow_rate_sphere_r','gas/density','gas/temperature','gas/metallicity',
#                                     'gas/magnetic_energy_specific','gas/cr_energy_specific'],
#                                     r=(0.2,'rvir'),dr=(0.01,'rvir'),weight=['gas/absmomentum_sphere_r','star/mass'],
#                                     remove_subs=True)
# proj.save_FITS('NUT_00035_mollweide.fits')
# plot_single_galaxy_healpix('NUT_00035_mollweide.fits',['gas/massflow_rate_sphere_r','gas/density','gas/temperature','gas/metallicity','gas/magnetic_energy_specific','gas/cr_energy_specific'])
