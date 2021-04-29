import matplotlib.pyplot as plt
import seaborn as sns
import seaborn as sns
sns.set(style="white")
plotting_dictionary = dict(
    density = {'cmap':sns.color_palette("mako", as_cmap=True),
                'text_over':'white',
                'label':r'$\rho$ [g/cm$^{-3}$]',
                'units':'g*cm**-3',
                'vmin':3e-27,
                'vmax':7e-22
    },
    temperature = {'cmap':sns.color_palette("rocket", as_cmap=True),
                    'text_over':'black',
                    'label':r'$T$ [K]',
                    'units':'K',
                    'vmin':5e+2,
                    'vmax':1e+6
    },
    metallicity = {'cmap':sns.diverging_palette(145, 300, s=60, as_cmap=True),
                    'text_over':'black',
                    'label':r'$Z$ [$Z_{\odot}$]',
                    'units':'dimensionless',
                    'vmin':5e-5,
                    'vmax':8e-3
    },
    v_sphere_r = {'cmap':sns.color_palette("vlag", as_cmap=True),
                    'text_over':'black',
                    'label':r'$v_r$ [km/s]',
                    'units':'km*s**-1'
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