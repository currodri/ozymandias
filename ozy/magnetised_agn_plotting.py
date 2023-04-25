import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import ozy
from ozy.plot_settings import plotting_dictionary
from ozy.projections import do_projection
from ozy.profiles import compute_profile
import matplotlib.font_manager as fm
from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def agn_plots(group, ozy_file, path=None, use_defaults=True, quantities=None, quantity_dicts=None, presentable=False, global_pov='z', system_type='AGN', smooth=False):
    """
    Module for quickly plotting all necessary observables for magnetised agn experiments
    NOTES: 
    -> Variables are computed and plotted for each window in turn, this is designed to be quick and dirty
    -> These are not designed to be publication/ready, but presentable will make them look nice :)
    """
    
    width = 5 # standard for MNRAS
    if presentable:
        plt.rcParams['text.usetex'] = True
        plt.rcParams['lines.linewidth'] = 1
        font = {'family' : 'sans', 'weight' : 'normal', 'size'   : 8}
        plt.rcParams.update({"font.family": "serif", "pgf.rcfonts": False, "axes.unicode_minus": False})
        plt.rc('font', **font)
        plt.rcParams['axes.edgecolor'] = '0.0'  

    
    fig, ax = plt.subplots(nrows=2, ncols=3, gridspec_kw={'width_ratios': [1, 1, 1]})
    plt.gcf().set_size_inches(3*width, 2*width)

    obj = group.obj
    if use_defaults:

        quantity_dicts = {
            '0/0': {
                'type': 'projection',
                'pov': f'{global_pov}',
                'quantity': ['gas/density'],
                'weight': ['gas/counts'],
                'default_plot_dict': True,
                'colormap': None,
                'colorbar': True,
                'snapshot': True,
                'axis_labels': True,  
                'logscale': True,             
            },
            '1/0': {
                'type': 'projection',
                'pov': f'{global_pov}',
                'quantity': ['gas/temperature'],
                'weight': ['gas/counts'],
                'default_plot_dict': True,
                'colormap': None,
                'colorbar': True,
                'snapshot': True,
                'axis_labels': True,  
                'logscale': True,             
            },
            '0/1': {
                'type': 'projection',
                'pov': f'{global_pov}',
                'quantity': ['gas/magnetic_energy'],
                'weight': ['gas/counts'],
                'default_plot_dict': True,
                'colormap': None,
                'colorbar': True,
                'snapshot': True,
                'axis_labels': True,  
                'logscale': True,
                'text_color': 'w',             
            },
            '1/1': {
                'type': 'projection',
                'pov': f'{global_pov}',
                'quantity': ['gas/magnetic_divergence'],
                'weight': ['gas/counts'],
                'default_plot_dict': True,
                'colormap': None,
                'colorbar': True,
                'snapshot': True,
                'axis_labels': True,  
                'logscale': True,             
            },
            '0/2': {
                'type': 'projection',
                'pov': f'{global_pov}',
                'quantity': ['gas/beta'],
                'weight': ['gas/counts'],
                'default_plot_dict': True,
                'colormap': None,
                'colorbar': True,
                'snapshot': True,
                'axis_labels': True,  
                'logscale': True,
                'text_color': 'w',             
            },
            '1/2': {
                'type': 'profile',
                'pov': f'{global_pov}',
                'x_var': 'r_sphere',
                'quantity': [['gas/massflow_rate_sphere_r']],
                'weight': ['gas/counts'],
                'axis_labels': True,  
                'y_log': True,
                'x_log': False,  
                'y_abs': True,        
            },
        }

    for key in quantity_dicts:
        #Computing quantities
        if quantity_dicts[key]['type'] == 'projection':
            print(quantity_dicts[key]['quantity'][0])
            #proj = do_projection(group,quantity_dicts[key]['quantity'], weight=quantity_dicts[key]['weight'], pov=quantity_dicts[key]['pov'])
            proj = do_projection(group,quantity_dicts[key]['quantity'], weight=quantity_dicts[key]['weight'], pov=quantity_dicts[key]['pov'], rmax=(2, 'kpc'))

            width_x =  0.675*480
            width_y =  0.675*480
            ax_key1 = int(key.split('/')[0]) # which row
            ax_key2 = int(key.split('/')[1]) # which column
            field = quantity_dicts[key]['quantity'][0]

            ex = [-0.5*width_x,0.5*width_x,-0.5*width_y,0.5*width_y]
            ax[ax_key1, ax_key2].set_xlim([-0.5*width_x,0.5*width_y])
            ax[ax_key1, ax_key2].set_ylim([-0.5*width_x,0.5*width_y])
            ax[ax_key1, ax_key2].axes.xaxis.set_visible(False)
            ax[ax_key1, ax_key2].axes.yaxis.set_visible(False)
            ax[ax_key1, ax_key2].axis('off')

            if smooth:
                # We smooth with a Gaussian kernel with x_stddev=1 (and y_stddev=1)
                # It is a 9x9 array
                from astropy.convolution import Gaussian2DKernel, convolve
                kernel = Gaussian2DKernel(x_stddev=0.8)
                # astropy's convolution replaces the NaN pixels with a kernel-weighted
                # interpolation from their neighbors
                data = convolve((proj.data_maps[0][0][0]).T, kernel)
                # sigma=3
                # cImage = sp.ndimage.filters.gaussian_filter(hdul[h].data.T, sigma, mode='constant')
            else:
                data = (proj.data_maps[0][0][0]).T


            
            if quantity_dicts[key]['colormap'] != None:
                colormap = quantity_dicts[key]['colormap']
                text_color = quantity_dicts[key]['text_color']
            else:
                plotting_def = plotting_dictionary[field.split('/')[1]]
                colormap = plotting_def['cmap']
                text_color = plotting_def['text_over']
                label = plotting_def['label']

            full_varname = field.split('/')[1]
                           

            if quantity_dicts[key]['logscale']:
                data = np.log10(data)
            
            plot = ax[ax_key1, ax_key2].imshow(data, cmap=colormap, origin='lower', interpolation='nearest', extent=ex)

            if quantity_dicts[key]['colorbar']:
                cbaxes = inset_axes(ax[ax_key1, ax_key2], width="80%", height="5%", loc='lower center')
                cbar = fig.colorbar(plot, cax=cbaxes, orientation='horizontal')
                if quantity_dicts[key]['logscale']:
                    #cbar.set_label(plotting_def['label_log'],color=plotting_def['text_over'],fontsize=10,labelpad=-20, y=0.85,weight='bold')
                    plot_label = plotting_def['label_log']
                else:
                    #cbar.set_label(plotting_def['label'],color=plotting_def['text_over'],fontsize=10,labelpad=-20, y=1.25)
                    plot_label = plotting_def['label']
                cbar.ax.xaxis.label.set_font_properties(fm.FontProperties(weight='bold',size=5))
                cbar.ax.tick_params(axis='x', pad=-40,labelsize=10,labelcolor=plotting_def['text_over'])
                cbar.ax.tick_params(length=0,width=0)

            fontprops = fm.FontProperties(size=20,weight='bold')

            
            ax[ax_key1, ax_key2].text(0.03, 0.87, plot_label,#f'{full_varname}', # Print field name
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax[ax_key1, ax_key2].transAxes,
                                color=text_color, fontsize=15,fontweight='bold')

        elif  quantity_dicts[key]['type'] == 'profile':
            print(quantity_dicts[key]['quantity'][0][0])
            width_x =  0.675*480
            width_y =  0.675*480

            for q in quantity_dicts[key]['quantity']:
                prof = compute_profile(group, ozy_file, quantity_dicts[key]['x_var'], q, quantity_dicts[key]['weight'], region_type='sphere', recompute=True, save=False, nbins=100,rmin=(0., 'code_length'), rmax=(0.5, 'code_length'))

                new_x = 0.5*(prof.xdata[0][:-1]+prof.xdata[0][1:])
                
                ax_key1 = int(key.split('/')[0]) # which row
                ax_key2 = int(key.split('/')[1]) # which column
                field = q[0]


                if quantity_dicts[key]['y_abs']:
                    ax[ax_key1, ax_key2].plot(new_x.in_units('kpc'), np.abs(prof.ydata['hydro'][0][:,0,0]), c='k', label=f"{field.split('/')[1]}")
                else:
                    ax[ax_key1, ax_key2].plot(new_x.in_units('kpc'), prof.ydata['hydro'][0][:,0,0], c='k', label=f"{field.split('/')[1]}")


            if quantity_dicts[key]['y_log']:
                ax[ax_key1, ax_key2].set_yscale('log')

            if quantity_dicts[key]['x_log']:
                ax[ax_key1, ax_key2].set_xscale('log')

            if len(quantity_dicts[key]['quantity']) > 1:
                ax[ax_key1, ax_key2].legend(frameon=False)
            else: 
                ax[ax_key1, ax_key2].set_ylabel(f"{quantity_dicts[key]['quantity'][0][0]}", rotation=270)

            ax[ax_key1, ax_key2].set_xlabel(f"{quantity_dicts[key]['x_var']}")

            ax[ax_key1, ax_key2].yaxis.set_label_position("right")
            ax[ax_key1, ax_key2].yaxis.tick_right()

        





    filename =  f'{system_type}_{obj.simulation.fullpath[-5:]}_{global_pov}'
    if path != None:
        filename = path + filename

    filename += '.png'

    plt.subplots_adjust(wspace=0., hspace=0.)
    plt.savefig(filename, format='png', dpi=330)
