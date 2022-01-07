import os
import numpy as np
import ozy
from ozy.projections import do_projection,plot_single_galaxy_projection,do_healpix_projection,plot_single_galaxy_healpix


obj = ozy.load('test_00035.hdf5')
virial_mass = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(virial_mass)
gal = obj.galaxies[progind]

# proj = do_projection(gal,['gas/density','gas/temperature','gas/metallicity','gas/magnetic_magnitude','star/sdensity','dm/sdensity'],window=(10,'kpc'),pov='faceon')
# proj.save_FITS('NUT_00010_faceon.fits')
# plot_single_galaxy_projection('NUT_00010_faceon.fits',['gas/density','gas/temperature','gas/metallicity','gas/magnetic_magnitude','star/sdensity','dm/sdensity'])


proj = do_healpix_projection(gal,['gas/v_sphere_r','gas/density','gas/temperature','gas/metallicity'],r=(0.2,'rvir'),dr=(0.01,'rvir'))
proj.save_FITS('NUT_00035_mollweide.fits')
plot_single_galaxy_healpix('NUT_00035_mollweide.fits',['gas/v_sphere_r','gas/density','gas/temperature','gas/metallicity'])