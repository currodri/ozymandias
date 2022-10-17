import numpy as np
import ozy
from ozy.projections import do_projection,plot_single_galaxy_projection
from ozy.outflow_inflow import compute_flows

obj = ozy.load('test_00035.hdf5')
masses = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(masses)
gal = obj.galaxies[progind]

# Check by eye whether the projections for the flows make sense
proj = do_projection(gal,['gas/density','gas/temperature','gas/v_sphere_r','gas/metallicity'],window=(0.5,'rvir'),pov='edgeon',
                    filter_conds=['none','v_sphere_r/>/15/km*s**-1','v_sphere_r/<=/-20/km*s**-1'],
                    filter_name=['all','outflow','inflow'],weight=['gas/volume','star/cumulative'])
proj.save_FITS('NUT_00035_edgeon.fits')
plot_single_galaxy_projection('NUT_00035_edgeon.fits',['gas/density','gas/temperature','gas/v_sphere_r','gas/metallicity'],
                                filter_name='all')
plot_single_galaxy_projection('NUT_00035_edgeon.fits',['gas/density','gas/temperature','gas/v_sphere_r','gas/metallicity'],
                                filter_name='outflow')
plot_single_galaxy_projection('NUT_00035_edgeon.fits',['gas/density','gas/temperature','gas/v_sphere_r','gas/metallicity'],
                                filter_name='inflow')

proj = do_projection(gal,['gas/density','gas/temperature','gas/v_sphere_r','gas/metallicity'],window=(0.5,'rvir'),pov='faceon',
                    filter_conds=['none','v_sphere_r/>/15/km*s**-1','v_sphere_r/<=/-20/km*s**-1'],
                    filter_name=['all','outflow','inflow'],weight=['gas/volume','star/cumulative'])
proj.save_FITS('NUT_00035_faceon.fits')
plot_single_galaxy_projection('NUT_00035_faceon.fits',['gas/density','gas/temperature','gas/v_sphere_r','gas/metallicity'],
                                filter_name='all')
plot_single_galaxy_projection('NUT_00035_faceon.fits',['gas/density','gas/temperature','gas/v_sphere_r','gas/metallicity'],
                                filter_name='outflow')
plot_single_galaxy_projection('NUT_00035_faceon.fits',['gas/density','gas/temperature','gas/v_sphere_r','gas/metallicity'],
                                filter_name='inflow')



# Compute flows

# compute_flows(gal,'test_00035.hdf5','outflow',rmin=(0.19,'rvir'),rmax=(0.21,'rvir'),save=True,recompute=True)
# compute_flows(gal,'test_00035.hdf5','inflow',rmin=(0.49,'rvir'),rmax=(0.51,'rvir'),save=True,recompute=True)