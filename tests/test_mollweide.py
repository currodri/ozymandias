import os
import numpy as np
import ozy
from ozy.projections import do_healpix_projection,plot_single_galaxy_healpix


obj = ozy.load('test_00010.hdf5')
virial_mass = [i.virial_quantities['mass'] for i in obj.galaxies]
progind = np.argmax(virial_mass)
gal = obj.galaxies[progind]


proj = do_healpix_projection(gal,['gas/massflow_rate_sphere_r','gas/density','gas/temperature','gas/metallicity',
                                    'gas/magnetic_energy_specific','gas/cr_energy_specific'],
                                    r=(1.0,'rvir'),dr=(0.7,'rvir'),weight=['gas/density','star/mass'],
                                    filter_conds=[[f'theta_sphere/<=/{np.pi/2}/radian'],
                                                  ['phi_sphere/>=/0.0/radian',
                                                   f'phi_sphere/</{np.pi}/radian']],
                                    filter_name=['top','posx'])
proj.save_FITS('NUT_00010_mollweide.fits')
plot_single_galaxy_healpix('NUT_00010_mollweide.fits',
                           ['gas/massflow_rate_sphere_r',
                            'gas/density','gas/temperature',
                            'gas/metallicity','gas/magnetic_energy_specific',
                            'gas/cr_energy_specific'],
                           type_scale='_galaxy',
                           filter_name='top')
plot_single_galaxy_healpix('NUT_00010_mollweide.fits',
                           ['gas/massflow_rate_sphere_r',
                            'gas/density','gas/temperature',
                            'gas/metallicity','gas/magnetic_energy_specific',
                            'gas/cr_energy_specific'],
                           type_scale='_galaxy',
                           filter_name='posx')