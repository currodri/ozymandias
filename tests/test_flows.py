import numpy as np
import ozy
from ozy.projections import do_projection,plot_single_galaxy_projection
from ozy.outflow_inflow import compute_flows

obj = ozy.load('test_00035.hdf5')
masses = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(masses)
gal = obj.galaxies[progind]

# Check by eye whether the projections for the flows make sense
# proj = do_projection(gal,['gas/density','gas/temperature','gas/v_sphere_r','gas/gradscale_crprsphere',
#                             'gas/grav_therpfrsphere','gas/grav_crpfrsphere'],window=(0.25,'rvir'),pov='edgeon',
#                             filter_conds=[['v_sphere_r/>/0/km*s**-1','temperature/>/1e5/K'],
#                                         ['v_sphere_r/>/0/km*s**-1','temperature/</1e5/K','temperature/>/9e3/K'],
#                                         ['v_sphere_r/>/0/km*s**-1','temperature/</9e3/K','temperature/>/1e3/K'],
#                                         ['v_sphere_r/>/0/km*s**-1','temperature/</1e3/K']],
#                             filter_name=['hot','warm_ionised','warm_neutral','cold'],weight=['gas/massflux_rate_sphere_r','star/cumulative'])
# proj.save_FITS('NUT_00035_edgeon.fits')
plot_single_galaxy_projection('NUT_00035_edgeon.fits',['gas/density','gas/temperature','gas/v_sphere_r','gas/gradscale_crprsphere',
                                'gas/grav_therpfrsphere','gas/grav_crpfrsphere'],
                                filter_name='hot')
plot_single_galaxy_projection('NUT_00035_edgeon.fits',['gas/density','gas/temperature','gas/v_sphere_r','gas/gradscale_crprsphere',
                                'gas/grav_therpfrsphere','gas/grav_crpfrsphere'],
                                filter_name='cold')
plot_single_galaxy_projection('NUT_00035_edgeon.fits',['gas/density','gas/temperature','gas/v_sphere_r','gas/gradscale_crprsphere',
                                'gas/grav_therpfrsphere','gas/grav_crpfrsphere'],
                                filter_name='warm_ionised')
plot_single_galaxy_projection('NUT_00035_edgeon.fits',['gas/density','gas/temperature','gas/v_sphere_r','gas/gradscale_crprsphere',
                                'gas/grav_therpfrsphere','gas/grav_crpfrsphere'],
                                filter_name='warm_neutral')

# proj = do_projection(gal,['gas/density','gas/temperature','gas/v_sphere_r','gas/metallicity'],window=(0.25,'rvir'),pov='faceon',
#                     filter_conds=['none','v_sphere_r/>/15/km*s**-1','v_sphere_r/<=/-20/km*s**-1'],
#                     filter_name=['all','outflow','inflow'],weight=['gas/density','star/cumulative'])
# proj.save_FITS('NUT_00035_faceon.fits')
# plot_single_galaxy_projection('NUT_00035_faceon.fits',['gas/density','gas/temperature','gas/v_sphere_r','gas/metallicity'],
#                                 filter_name='all')
# plot_single_galaxy_projection('NUT_00035_faceon.fits',['gas/density','gas/temperature','gas/v_sphere_r','gas/metallicity'],
#                                 filter_name='outflow')
# plot_single_galaxy_projection('NUT_00035_faceon.fits',['gas/density','gas/temperature','gas/v_sphere_r','gas/metallicity'],
#                                 filter_name='inflow')



# Compute flows

# flow = compute_flows(gal,'test_00035.hdf5','outflow',rmin=(0.195,'rvir'),rmax=(0.205,'rvir'),save=True,recompute=True,remove_subs=False)
#print(flow.data)
# compute_flows(gal,'test_00035.hdf5','inflow',rmin=(0.49,'rvir'),rmax=(0.51,'rvir'),save=True,recompute=True)