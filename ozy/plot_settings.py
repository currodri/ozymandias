import matplotlib.pyplot as plt
import swiftascmaps
import seaborn as sns
sns.set(style="white")
plotting_dictionary = dict(
    density = {'cmap': sns.color_palette("mako", as_cmap=True),
                'text_over':'white',
                'label':r'$\rho$ [g/cm$^{-3}$]',
                'label_log':r'$\log\left(\frac{\rho}{{\rm g cm}^{-3}}\right)$',
                'units':'g*cm**-3',
                'vmin':2e-30,
                'vmax':7e-22,
                'vmin_galaxy':2e-28,
                'vmax_galaxy':5e-22
    },
    mass = {'cmap':sns.color_palette("mako", as_cmap=True),
            'text_over':'white',
            'label':r'Mass [M$_{\odot}$]',
            'label_log':r'$\log$(\frac{Mass}{{\rmM}$_{\odot}$})',
            'units':'Msun',
            'vmin':1,
            'vmax':1e+8
    },
    temperature = {'cmap':sns.color_palette("rocket", as_cmap=True),
                    'text_over':'white',
                    'label':r'$T$ [K]',
                    'label_log':r'$\log\left(\frac{T}{{\rm K}}\right)$',
                    'units':'K',
                    'vmin':2e+4,
                    'vmax':8e+6,
                    'vmin_galaxy':5e+2,
                    'vmax_galaxy':8e+5
    },
    thermal_energy = {'cmap':sns.color_palette("rocket", as_cmap=True),
                            'text_over':'white',
                            'label':r'$E_{\rm ther}$ [erg]',
                            'label_log':r'$\log(\frac{E_{\rm ther}}{{\rm erg}})$',
                            'units':'erg',
                            'vmin':5e+1,
                            'vmax':2e+6,
                            'vmin_galaxy':5e+3,
                            'vmax_galaxy':8e+6
    },
    thermal_energy_specific = {'cmap':sns.color_palette("rocket", as_cmap=True),
                                'text_over':'white',
                                'label':r'$\epsilon_{\rm ther}$ [erg/g]',
                                'label_log':r'$\log(\frac{\epsilon_{\rm ther}}{{\rm erg/g}})$',
                                'units':'erg*g**-1',
                                'vmin':5e+1,
                                'vmax':2e+6,
                                'vmin_galaxy':5e+3,
                                'vmax_galaxy':8e+6
    },
    entropy_specific = {'cmap':'plasma',
                        'text_over':'white',
                        'label':r'$s$ [erg/(g K)]',
                        'label_log':r'$\log(\frac{s}{{\rm erg/(g K)}})$',
                        'units':'erg*K**-1*g**-1',
                        'vmin':8e+3,
                        'vmax':2e+5,
                        'vmin_galaxy':8e+3,
                        'vmax_galaxy':2e+5
    },
    metallicity = {'cmap':'swift.nineteen_eighty_nine',
                    'text_over':'black',
                    'label':r'$Z$ [$Z_{\odot}$]',
                    'label_log':r'$\log\left(Z/Z_{\odot}\right)$',
                    'units':'dimensionless',
                    'vmin':5e-3,
                    'vmax':1.5,
                    'vmin_galaxy':5e-3,
                    'vmax_galaxy':1.5
    },
    magnetic_energy = {'cmap':sns.cubehelix_palette(reverse=True,as_cmap=True),
                                'text_over':'white',
                                'label':r'$E_{\rm mag}$ [erg]',
                                'label_log':r'$\log(\frac{E_{\rm mag}}{{\rm erg}})$',
                                'units':'erg',
                                'vmin':5e+10,
                                'vmax':3e+12

    },
    magnetic_energy_specific = {'cmap':sns.cubehelix_palette(reverse=True,as_cmap=True),
                                'text_over':'white',
                                'label':r'$\epsilon_{\rm mag}$ [erg/g]',
                                'units':'erg*g**-1',
                                'vmin':5e+8,
                                'vmax':3e+11,
                                'vmin_galaxy':5e+8,
                                'vmax_galaxy':3e+11

    },
    magnetic_energy_density = {'cmap':sns.cubehelix_palette(reverse=True,as_cmap=True),
                                'text_over':'white',
                                'label':r'$\varepsilon_{\rm mag}$ [erg/cm$^{3}$]',
                                'label_log':r'$\log(\frac{\varepsilon_{\rm mag}}{{\rm erg/cm}^{3}})$',
                                'units':'erg*cm**-3',
                                'vmin':4e-17,
                                'vmax':5e-11,
                                'vmin_galaxy':4e-17,
                                'vmax_galaxy':5e-11

    },
    magnetic_magnitude = {'cmap':'inferno',
                                'text_over':'white',
                                'label':r'$B$ [$G$]',
                                'label_log':r'$\log(\frac{B}{{\rm G}})$',
                                'units':'G',
                                'vmin':4e-17,
                                'vmax':5e-11,
                                'vmin_galaxy':4e-10,
                                'vmax_galaxy':5e-5

    },
    alfven_speed = {'cmap':'swift.red',
                                'text_over':'white',
                                'label':r'$v_A$ [km s$^{-1}$]',
                                'label_log':r'$\log(\frac{v_A}{{\rm km s}^{-1}})$',
                                'units':'km/s',
                                'vmin':3e-2,
                                'vmax':5e0,
                                'vmin_galaxy':3e-2,
                                'vmax_galaxy':5e0

    },
    total_crs_energy = {'cmap':sns.cubehelix_palette(start=2, rot=0, dark=0, light=.95, reverse=True, as_cmap=True),
                        'text_over':'white',
                        'label':r'$E_{\rm CR}$ [erg]',
                        'label_log':r'$\log(\frac{\varepsilon_{\rm CR}}{{\rm erg/cm}^{3}})$',
                        'label_log':r'$\log(\frac{E_{\rm CR}}{{\rm erg}})$',
                        'units':'erg',
                        'vmin':4e-14,
                        'vmax':8e-11,
                        'vmin_galaxy':4e-14,
                        'vmax_galaxy':8e-11

    },
    cr_energy_density = {'cmap':sns.cubehelix_palette(start=2, rot=0, dark=0, light=.95, reverse=False, as_cmap=True),
                        'text_over':'white',
                        'label':r'$\varepsilon_{\rm CR}$ [erg/cm$^{3}$]',
                        'label_log':r'$\log(\frac{\varepsilon_{\rm CR}}{{\rm erg/cm}^{3}})$',
                        'units':'erg*cm**-3',
                        'vmin':7e-18,
                        'vmax':8e-12,
                        'vmin_galaxy':8e-14,
                        'vmax_galaxy':8e-11

    },
    cr_energy_specific = {'cmap':sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True), #sns.cubehelix_palette(start=2, rot=0, dark=0, light=.95, reverse=True, as_cmap=True),
                                'text_over':'white',
                                'label':r'$\epsilon_{\rm CR}$ [erg/g]',
                                'label_log':r'$\log(\frac{\epsilon_{\rm CR}}{\rm erg/g})$',
                                'units':'erg*g**-1',
                                'vmin':5e+12,
                                'vmax':7e+14,
                                'vmin_galaxy':5e+10,
                                'vmax_galaxy':7e+15

    },
    xHII = {'cmap':'bone',
            'text_over':'white',
            'label':r'$X_{\rm HII}$',
            'label_log':r'$\log(X_{\rm HII})$',
            'units':'dimensionless',
            'vmin':5e-4,
            'vmax':0.8,
            'vmin_galaxy':5e-4,
            'vmax_galaxy':0.8

    },
    xHeII = {'cmap':'swift.red',
            'text_over':'white',
            'label':r'$X_{\rm HeII}$',
            'label_log':r'$\log(X_{\rm HeII})$',
            'units':'dimensionless',
            'vmin':5e-4,
            'vmax':0.8,
            'vmin_galaxy':5e-4,
            'vmax_galaxy':0.8

    },
    xHeIII = {'cmap':'swift.red',
                'text_over':'white',
                'label':r'$X_{\rm HeIII}$',
                'label_log':r'$\log(X_{\rm HeIII})$',
                'units':'dimensionless',
                'vmin':5e-4,
                'vmax':0.8

    },
    v_sphere_r = {'cmap':sns.color_palette("vlag", as_cmap=True),
                    'text_over':'black',
                    'label':r'$v_r$ [km/s]',
                    'label_log':r'$\log(\frac{v_r}{\rm km/s})$',
                    'units':'km*s**-1',
                    'vmin':-90,
                    'vmax':+90
    },
    momentum_sphere_r = {'cmap':sns.color_palette("vlag", as_cmap=True),
                    'text_over':'black',
                    'label':r'$p_r$ [M$_{\odot}$ km/s]',
                    'units':'Msun*km*s**-1'
    },
    massflow_rate = {'cmap':sns.color_palette("vlag", as_cmap=True),
                                'text_over':'black',
                                'label':r'$dM/dt$ [M$_{\odot}$ yr$^{-1}$]',
                                'units':'Msun*yr**-1'
    },
    massflow_rate_sphere_r = {'cmap':sns.color_palette("vlag", as_cmap=True),
                                'text_over':'black',
                                'label':r'$dM/dt$ [M$_{\odot}$ yr$^{-1}$]',
                                'units':'Msun*yr**-1'
    },
    star_mass = {'cmap':'gray',
                'text_over':'white',
                'label':r'$M_{*}$ [M$_{\odot}$]',
                'label_log':r'$\log\left(\frac{M_{*}}{{\rmM}_{\odot}}\right)$',
                'units':'Msun',
                'vmin':5.0,
                'vmax':2e+7
    },
    star_sdensity = {'cmap':'gray',
                        'text_over':'white',
                        'label':r'$\Sigma_{*}$ [M$_{\odot}$ kpc$^{-2}$]',
                        'label_log':r'$\log\left(\frac{\Sigma_{*}}{{\rmM}_{\odot}{\rm kpc}^{-2}}\right)$',
                        'units':'Msun/(kpc**2)',
                        'vmin':3e+4,
                        'vmax':5e+9,
                        'vmin_galaxy':3e+4,
                        'vmax_galaxy':5e+9

    },
    dm_sdensity = {'cmap':'cividis',
                        'text_over':'white',
                        'label':r'$\Sigma_{\rm DM}$ [M$_{\odot}$ kpc$^{-2}$]',
                        'label_log':r'$\log\left(\frac{\Sigma_{\rm DM}}{{\rmM}_{\odot}{\rm kpc}^{-2}}\right)$',
                        'units':'Msun/(kpc**2)',
                        # 'vmin':4e+5,
                        # 'vmax':8e+8
                        'vmin':4e+7,
                        'vmax':3e+9

    },
    dm_mass = {'cmap':'cividis',
                'text_over':'white',
                'label':r'$M_{DM}$ [M$_{\odot}$]',
                'label_log':r'$\log\left(\frac{M_{\rm DM}}{{\rmM}_{\odot}}\right)$',
                'units':'Msun',
                'vmin':3.0e+3,
                'vmax':7e+6
    },
    star_metallicity = {'cmap':sns.color_palette("dark:salmon", as_cmap=True),
                        'text_over':'white',
                        'label':r'$Z$ [$Z_{\odot}$]',
                        'units':'dimensionless',
                        'vmin':1e-1,
                        'vmax':None
    },
    star_age = {'cmap':'BuPu',
                'text_over':'white',
                'label':r'Age [Gyr]',
                'units':'Gyr',
                'vmin':None,
                'vmax':None
    },
    star_sfr_surface_100 = {'cmap':'magma',
                            'text_over':'white',
                            'label':r'$\Sigma_{\rm SFR}$ [M$_{\odot}$ yr$^{-1}$ kpc$^{-2}$]',
                            'label_log':r'$\log\left(\frac{\Sigma_{\rm SFR}}{{\rmM}_{\odot}{\rm yr}^{-1}{\rm kpc}^{-2}}\right)$',
                            'units':'Msun/(yr*kpc**2)',
                            'vmin':3e-4,
                            'vmax':90,
                            'vmin_galaxy':3e-4,
                            'vmax_galaxy':90
    }
)

circle_dictionary = dict(
        rvir_galaxy = {'edgecolor': 'y',
                        'linestyle':'--'
        },
        rvir_satellite = {'edgecolor': 'w',
                          'linestyle':'--'
        },
        tidal_BT87_simple = {'edgecolor': 'r',
                          'linestyle':':'
        },
        tidal_BT87_centrifugal = {'edgecolor': 'b',
                          'linestyle':':'
        }
)