import matplotlib.pyplot as plt
import seaborn as sns
import seaborn as sns
sns.set(style="white")
plotting_dictionary = dict(
    density = {'cmap':sns.color_palette("mako", as_cmap=True),
                'text_over':'white',
                'label':r'$\rho$ [g/cm$^{-3}$]',
                'label_log':r'$\log(\rho$ [g/cm$^{-3}$])',
                'units':'g*cm**-3',
                'vmin':3e-27,
                'vmax':3e-21
    },
    mass = {'cmap':sns.color_palette("mako", as_cmap=True),
            'text_over':'white',
            'label':r'Mass [M$_{\odot}$]',
            'label_log':r'$\log$(Mass) [M$_{\odot}$]',
            'units':'Msun',
            'vmin':1,
            'vmax':1e+8
    },
    temperature = {'cmap':sns.color_palette("rocket", as_cmap=True),
                    'text_over':'black',
                    'label':r'$T$ [K]',
                    'label_log':r'$\log(T$ [K])',
                    'units':'K',
                    'vmin':5e+1,
                    'vmax':2e+6
    },
    metallicity = {'cmap':sns.diverging_palette(145, 300, s=60, as_cmap=True),
                    'text_over':'black',
                    'label':r'$Z$ [$Z_{\odot}$]',
                    'units':'dimensionless',
                    'vmin':6e-4,
                    'vmax':9e-1
    },
    magnetic_energy_specific = {'cmap':sns.cubehelix_palette(reverse=True,as_cmap=True),
                                'text_over':'white',
                                'label':r'$\epsilon_{\rm mag}$ [erg/g]',
                                'units':'erg*g**-1',
                                'vmin':1e+10,
                                'vmax':1e+12

    },
    magnetic_energy_density = {'cmap':sns.cubehelix_palette(reverse=True,as_cmap=True),
                                'text_over':'white',
                                'label':r'$\varepsilon_{\rm mag}$ [erg/cm$^{3}$]',
                                'units':'erg*g**-1',
                                'vmin':4e-19,
                                'vmax':1e-10

    },
    v_sphere_r = {'cmap':sns.color_palette("vlag", as_cmap=True),
                    'text_over':'black',
                    'label':r'$v_r$ [km/s]',
                    'units':'km*s**-1',
                    'vmin':-120,
                    'vmax':+120
    },
    momentum_sphere_r = {'cmap':sns.color_palette("vlag", as_cmap=True),
                    'text_over':'black',
                    'label':r'$p_r$ [M$_{\odot}$ km/s]',
                    'units':'Msun*km*s**-1'
    },
    star_mass = {'cmap':'gray',
                'text_over':'white',
                'label':r'$M_{*}$ [M$_{\odot}$]',
                'units':'Msun',
                'vmin':1e+5,
                'vmax':1e+8
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
    }
)