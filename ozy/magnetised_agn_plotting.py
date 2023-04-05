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


def agn_plots(group, path=None, use_defaults=True, quantities=None, quantity_dicts=None, presentable=False):
    """
    Module for quickly plotting all necessary observables for magnetised agn experiments
    NOTES: 
    -> Variables are computed and plotted for each window in turn, this is designed to be quick and dirty
    -> These are not designed to be publication/ready, but presentable will make them look nice :)
    """
    
    width = 3.31 # standard for MNRAS
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
        global_pov = 'y' #TODO: Make this a user-parameter

        quantity_dicts = {
            '0/0': {
                'type': 'projection',
                'pov': f'{global_pov}',
                'quantity': ['gas/density'],
                'weight': ['gas/counts'],
                'default_plot_dict': True,
                'colormap': None,
                'colorbar': False,
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
                'colorbar': False,
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
                'colormap': 'plasma',
                'colorbar': False,
                'snapshot': True,
                'axis_labels': True,  
                'logscale': True,
                'text_color': 'w',             
            },
            '1/1': {
                'type': 'projection',
                'pov': f'{global_pov}',
                'quantity': ['gas/alfven_speed'],
                'weight': ['gas/counts'],
                'default_plot_dict': True,
                'colormap': None,
                'colorbar': False,
                'snapshot': True,
                'axis_labels': True,  
                'logscale': True,             
            },
            '0/2': {
                'type': 'projection',
                'pov': f'{global_pov}',
                'quantity': ['gas/momentum_sphere_r'],
                'weight': ['gas/counts'],
                'default_plot_dict': True,
                'colormap': 'icefire',
                'colorbar': False,
                'snapshot': True,
                'axis_labels': True,  
                'logscale': False,
                'text_color': 'w',             
            },
            '1/2': {
                'type': 'projection',
                'pov': f'{global_pov}',
                'quantity': ['gas/thermal_pressure'],
                'weight': ['gas/counts'],
                'default_plot_dict': True,
                'colormap': 'bone',
                'colorbar': False,
                'snapshot': True,
                'axis_labels': True,  
                'logscale': True,
                'text_color': 'w',          
            },
        }

    for key in quantity_dicts:
        #Computing quantities
        if quantity_dicts[key]['type'] == 'projection':
            print(quantity_dicts[key]['quantity'][0])
            proj = do_projection(group,quantity_dicts[key]['quantity'], weight=quantity_dicts[key]['weight'], pov=quantity_dicts[key]['pov'])

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
                    cbar.set_label(plotting_def['label_log'],color=plotting_def['text_over'],fontsize=10,labelpad=-60, y=0.85,weight='bold')
                else:
                    cbar.set_label(plotting_def['label'],color=plotting_def['text_over'],fontsize=10,labelpad=-10, y=1.25)
                cbar.ax.xaxis.label.set_font_properties(fm.FontProperties(weight='bold',size=5))
                cbar.ax.tick_params(axis='x', pad=-16, labelsize=13,labelcolor=plotting_def['text_over'])
                cbar.ax.tick_params(length=0,width=0)

            fontprops = fm.FontProperties(size=20,weight='bold')

            
            ax[ax_key1, ax_key2].text(0.03, 0.87, f'{full_varname}', # Print field name
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax[ax_key1, ax_key2].transAxes,
                                color=text_color, fontsize=10,fontweight='bold')

        elif  quantity_dicts[key]['type'] == 'profile':
            print(quantity_dicts[key]['quantity'][0])
            #compute_profile()


    filename =  f'AGN_{obj.simulation.fullpath[-5:]}_{global_pov}'
    if path != None:
        filename = path + filename

    filename += '.png'

    plt.subplots_adjust(wspace=0.)
    plt.savefig(filename, format='png', dpi=330)