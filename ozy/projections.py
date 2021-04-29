import numpy as np
import h5py
import os
import ozy
from ozy.plot_settings import plotting_dictionary
from ozy.dict_variables import common_variables,grid_variables,particle_variables,get_code_units,basic_conv
import sys
sys.path.append('/mnt/extraspace/currodri/Codes/ozymandias/ozy/visualisation')
import re
from projections import obs_instruments
from projections import maps
from projections import vectors
from astropy.io import fits
from astropy.wcs import WCS


cartesian_basis = {'x':np.array([1.,0.,0.]),'y':np.array([0.,0.,1.]),'z':np.array([0.,0.,1.])}

class Projection(object):

    def __init__(self,group):
        self.group = group
        self.pov = 'faceon'
        self.los_axis = np.array([0.,0.,1.])
        self.up_vector = np.array([0.,1.,0.])
        self.center = np.array([0.5,0.5,0.5])
        self.width = np.array([1.,1.])
        self.distance = 0.5
        self.far_cut_depth = 0.5
        self.resolution = np.array([1024,1024],'i')
        self.vars = dict(gas = [],star = [], dm = [])
        self.weight = ['density','star/density']
        self.data_maps = []

    def save_FITS(self,name,unit_system=basic_conv):
        if os.path.exists(name):
            print('WARNING: Overwritting previous projection.')
            os.remove(name)
        from yt import YTArray, YTQuantity
        self.hdulist = fits.HDUList()
        first = True
        counter = 0
        for datatype, varlist in self.vars.items():
            fields = []

            if len(varlist)>0:
                imap = self.data_maps[counter]
                for f in varlist:
                    fields.append(datatype+'/'+f)
                for i,field in enumerate(fields):
                    code_units = get_code_units(field.split('/')[1])
                    temp_map = YTArray(imap[i],code_units,
                                        registry=self.group.obj.unit_registry)
                    first_unit = True
                    for u in code_units.split('*'):
                        if first_unit:
                            units = unit_system[u]
                        else:
                            units += '*'+unit_system[u]
                    if first:
                        hdu = fits.PrimaryHDU(np.array(temp_map.in_units(units)))
                        first = False
                    else:
                        hdu = fits.ImageHDU(np.array(temp_map.in_units(units)))
                    hdu.name = field
                    hdu.header["btype"] = field
                    hdu.header["bunit"] = re.sub('()', '', units)
                    hdu.header["redshift"] = self.group.obj.simulation.redshift
                    self.hdulist.append(hdu)
            counter += 1
        # Setup WCS in the coordinate systems of the camera
        w = WCS(header=self.hdulist[0].header,naxis=2)
        dx = YTQuantity(self.width[0],'code_length',
                        registry=self.group.obj.unit_registry).in_units(unit_system['code_length']).d
        dy = YTQuantity(self.width[1],'code_length',
                        registry=self.group.obj.unit_registry).in_units(unit_system['code_length']).d
        dx /= self.resolution[0]
        dy /= self.resolution[1]
        center = [0.0,0.0]
        cdelt = [dx,dy]
        w.wcs.crpix = 0.5*(np.array(self.resolution)+1)
        w.wcs.crval = center
        w.wcs.cdelt = cdelt
        w.wcs.ctype = ["linear"]*2
        w.wcs.cunit = [unit_system['code_length']]*2
        self.set_wcs(w)

        #TODO: Add further information such that the FITS file is a standalone.

        # Now, just save to disk
        self.hdulist.writeto(name)

    def set_wcs(self,wcs):
        """Set the WCS coordinate information for all images in a Projection class."""
        wcs.wcs.name = self.pov
        h = wcs.to_header()
        self.wcs = wcs
        for img in self.hdulist:
            for k,v in h.items():
                kk = k
                img.header[kk] = v

def do_projection(group,vars,weight=['gas/density','star/age'],map_max_size=1024,pov='faceon'):
    """Function which computes a 2D projection centered on an objected from an OZY file."""

    if not isinstance(weight,list):
        weight = [weight]
    if len(weight)>2:
        print('Only just one weight per data type (i.e. grid and particles)!')
        exit
    proj = Projection(group)
    proj.pov = pov
    
    for var in vars:
        var_type = var.split('/')[0]
        var_name = var.split('/')[1]
        if var_type == 'gas':
            if var_name in common_variables or var_name in grid_variables:
                proj.vars['gas'].append(var_name)
            else:
                raise KeyError('This gas variable is not supported. Please check!: %s', var)
        elif var_type == 'star':
            if var_name in common_variables or var_name in particle_variables:
                proj.vars['star'].append(var_name)
            else:
                raise KeyError('This star variable is not supported. Please check!')
        elif var_type == 'dm':
            if var_name in common_variables or var_name in particle_variables:
                proj.vars['dm'].append(var_name)
            else:
                raise KeyError('This DM variable is not supported. Please check!')

    for w in weight:
        weight_type = w.split('/')[0]
        weight_name = w.split('/')[1]
        if var_type == 'gas':
            if weight_name in common_variables or weight_name in grid_variables:
                proj.weight[0] = weight_name
            else:
                raise KeyError('This gas variable is not supported. Please check!: %s', var)
        elif var_type == 'star':
            if weight_name in common_variables or weight_name in particle_variables:
                proj.weight[1] = w
            else:
                raise KeyError('This star variable is not supported. Please check!')
        elif var_type == 'dm':
            if weight_name in common_variables or weight_name in particle_variables:
                proj.weight[1] = w
            else:
                raise KeyError('This DM variable is not supported. Please check!')

    # Setup camera details for the requested POV (Point of View)
    centre = vectors.vector()
    centre.x, centre.y, centre.z = group.position[0], group.position[1], group.position[2]
    bulk = vectors.vector()
    bulk.x, bulk.y, bulk.z = group.velocity[0], group.velocity[1], group.velocity[2]
    if proj.pov == 'faceon':
        axis = vectors.vector()
        norm_L = group.angular_mom['total'].d/np.linalg.norm(group.angular_mom['total'].d)
        up = cartesian_basis['y'] - np.dot(cartesian_basis['y'],norm_L)*norm_L
        up /= np.linalg.norm(up)
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        up_vector = vectors.vector()
        up_vector.x, up_vector.y, up_vector.z = up[0], up[1], up[2]
        rmax = 0.2*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        region_size = np.array([2*rmax,2*rmax],order='F')
        distance = rmax
        far_cut_depth = rmax
    elif proj.pov == 'edgeon':
        axis = vectors.vector()
        norm_L = group.angular_mom['total'].d/np.linalg.norm(group.angular_mom['total'].d)
        los = cartesian_basis['y'] - np.dot(cartesian_basis['y'],norm_L)*norm_L
        los /= np.linalg.norm(los)
        axis.x,axis.y,axis.z = los[0], los[1], los[2]
        up_vector = vectors.vector()
        up_vector.x,up_vector.y,up_vector.z = norm_L[0], norm_L[1], norm_L[2]
        rmax = 0.2*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        region_size = np.array([2*rmax,2*rmax],order='F')
        distance = rmax
        far_cut_depth = rmax
    else:
        print("This point of view is not supported!")
        print("Falling back to 'faceon' (default).")
        axis = vectors.vector()
        norm_L = group.angular_mom['total']/np.linalg.norm(group.angular_mom['total'])
        up = cartesian_basis['y'] - np.dot(cartesian_basis['y'],norm_L)
        up /= np.linalg.norm(up)
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        up_vector = vectors.vector()
        up_vector.x, up_vector.y, up_vector.z = up[0], up[1], up[2]
        rmax = 0.2*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        region_size = np.array([2*rmax,2*rmax],order='F')
        distance = rmax
        far_cut_depth = rmax
    
    cam = obs_instruments.init_camera(centre,axis,up_vector,region_size,distance,far_cut_depth,map_max_size-1)

    # Update projection details with the camera ones
    proj.width = region_size
    proj.centre = np.array([group.position[0], group.position[1], group.position[2]])
    proj.distance = distance
    proj.far_cut_depth = far_cut_depth
    proj.los_axis = norm_L
    proj.up_vector = np.array([cam.up_vector.x,cam.up_vector.y,cam.up_vector.z])

    # Save the final resolution that will be computed (such as for aspect_ratio>1)
    obs_instruments.get_map_size(cam,proj.resolution)

    # Create projection_handler Fortran derived type for the results of the hydro data projection
    hydro_handler = maps.projection_handler()
    hydro_handler.type = proj.pov
    hydro_handler.nvars = len(proj.vars['gas'])
    hydro_handler.weightvar = proj.weight[0]
    maps.allocate_projection_handler(hydro_handler)
    for i in range(0, len(proj.vars['gas'])):
        hydro_handler.varnames.T.view('S128')[i] = proj.vars['gas'][i].ljust(128)
    
    # COMPUTE HYDRO PROJECTION
    maps.projection_hydro(group.obj.simulation.fullpath,cam,bulk,hydro_handler)
    # TODO: Weird issue when the direct toto array is given.
    data = np.copy(hydro_handler.toto)
    proj.data_maps.append(data)

    # Create projection_handler Fortran derived type for the results of the hydro data projection
    parts_handler = maps.projection_handler()
    parts_handler.type = proj.pov
    parts_handler.nvars = len(proj.vars['dm'])+len(proj.vars['star'])
    parts_handler.weightvar = proj.weight[1]
    maps.allocate_projection_handler(parts_handler)
    for i in range(0, len(proj.vars['dm'])):
        tempstr = 'dm/'+proj.vars['dm'][i]
        parts_handler.varnames.T.view('S128')[i] = tempstr.ljust(128)
    for i in range(0, len(proj.vars['star'])):
        tempstr = 'star/'+proj.vars['star'][i]
        parts_handler.varnames.T.view('S128')[i] = tempstr.ljust(128)

    # COMPUTE PARTICLES PROJECTION
    maps.projection_parts(group.obj.simulation.fullpath,cam,bulk,parts_handler)
    # TODO: Weird issue when the direct toto array is given.
    data = np.copy(parts_handler.toto)
    proj.data_maps.append(data)

    return proj



def plot_single_galaxy_projection(proj_FITS,fields,logscale=True,scalebar=True,redshift=True):
    """This function uses the projection information in a FITS file following the 
        OZY format and plots it following the OZY standards."""
    
    # Make required imports
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from matplotlib.colors import LogNorm
    import seaborn as sns
    sns.set(style="white")
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    hfont = {'fontname':'Helvetica'}
    matplotlib.rc('text', usetex = True)
    matplotlib.rc('font', **{'family' : "serif"})
    params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
    matplotlib.rcParams.update(params)


    # First,check that FITS file actually exists
    if not os.path.exists(proj_FITS):
        raise ImportError('File not found. Please check!')
    
    # Load FITS file
    hdul = fits.open(proj_FITS)
    hdul_fields = [h.header['btype'] for h in hdul]

    # Check that the required fields for plotting are in this FITS
    for i,f in enumerate(fields):
        if f not in hdul_fields:
            print('The field %s is not included in this file. Ignoring...'%f)
            del f[i]
    if len(fields) == 0:
        print('Not a single field of the ones you provided are here... Check!')
        exit

    # Since everything is fine, we begin plotting…
    ncolumns = len(fields)
    fig = plt.figure(figsize=(11,5))
    grid = AxesGrid(fig, (111),
                    nrows_ncols = (1, ncolumns),
                    axes_pad = 0.0,
                    label_mode = "L",
                    share_all = True,
                    cbar_location="top",
                    cbar_mode="edge",
                    cbar_size="2%",
                    cbar_pad=0.0)
    width_x =  hdul[0].header['CDELT1']*hdul[0].header['NAXIS1']
    width_y =  hdul[0].header['CDELT2']*hdul[0].header['NAXIS2']
    ex = [-0.5*width_x,0.5*width_x,-0.5*width_y,0.5*width_y]

    stellar = False
    for i in range(0, len(fields)):
        h = [j for j in range(0,len(hdul)) if hdul[j].header['btype']==fields[i]][0]
        ax = grid[i].axes
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
        if fields[i].split('/')[0] == 'star' or fields[i].split('/')[0] == 'dm':
            plotting_def = plotting_dictionary[fields[i].split('/')[0]+'_'+fields[i].split('/')[1]]
            stellar = True
        else:
            plotting_def = plotting_dictionary[fields[i].split('/')[1]]
        if logscale:
            plot = ax.imshow(hdul[h].data.T, cmap=plotting_def['cmap'],
                            origin='upper',norm=LogNorm(vmin=plotting_def['vmin'],
                            vmax=plotting_def['vmax']),extent=ex,
                            interpolation='nearest')
        else:
            plot = ax.imshow(hdul[h].data.T, cmap=plotting_def['cmap'],
                            origin='upper',extent=ex,interpolation='nearest',
                            vmin=plotting_def['vmin'],vmax=plotting_def['vmax'])
        fontprops = fm.FontProperties(size=14)
        if scalebar:
            scalebar = AnchoredSizeBar(ax.transData,
                                        3, '3 kpc', 'lower left', 
                                        pad=0.1,
                                        color=plotting_def['text_over'],
                                        frameon=False,
                                        size_vertical=0.2,
                                        fontproperties=fontprops)
            ax.add_artist(scalebar)
        axcb = fig.colorbar(plot, cax = grid.cbar_axes[i], orientation='horizontal')
        axcb.set_label(plotting_def['label'], fontsize=16,labelpad=-50, y=0.85)
        axcb.ax.xaxis.set_ticks_position("top")
        if redshift:
            ax.text(0.7, 0.12, 'z = '+str(round(hdul[h].header['redshift'], 2)),
                    transform=ax.transAxes, fontsize=18,verticalalignment='top',
                    color=plotting_def['text_over'])

    fig.subplots_adjust(hspace=0, wspace=0,top = 0.95,bottom = 0.02,left = 0.02,right = 0.98)
    if stellar:
        fig.savefig(proj_FITS.split('.')[0]+'_stars.png',format='png',dpi=300)
    else:
        fig.savefig(proj_FITS.split('.')[0]+'.png',format='png',dpi=300)