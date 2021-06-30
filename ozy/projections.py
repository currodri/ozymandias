import numpy as np
import h5py
import os
from yt import YTArray
import ozy
from ozy.plot_settings import plotting_dictionary
from ozy.dict_variables import common_variables,grid_variables,particle_variables,get_code_units,basic_conv
import sys
sys.path.append('/mnt/extraspace/currodri/Codes/ozymandias/ozy/visualisation')
import re
import healpy as hp
from astropy.io import fits
from astropy.wcs import WCS
from projections import obs_instruments
from projections import maps
from projections import vectors
from projections import geometrical_regions

cartesian_basis = {'x':np.array([1.,0.,0.]),'y':np.array([0.,0.,1.]),'z':np.array([0.,0.,1.])}

target_header = fits.Header.fromstring("""
NAXIS   =                    2
NAXIS1  =                  480
NAXIS2  =                  240
CTYPE1  = 'RA---MOL'
CRPIX1  =                240.5
CRVAL1  =                180.0
CDELT1  =               -0.675
CUNIT1  = 'deg     '
CTYPE2  = 'DEC--MOL'
CRPIX2  =                120.5
CRVAL2  =                  0.0
CDELT2  =                0.675
CUNIT2  = 'deg     '
COORDSYS= 'icrs    '
""", sep='\n')

class Projection(object):

    def __init__(self,group):
        self.group = group
        self.pov = 'faceon'
        self.los_axis = np.array([0.,0.,1.])
        self.up_vector = np.array([0.,1.,0.])
        self.centre = np.array([0.5,0.5,0.5])
        self.width = np.array([1.,1.])
        self.distance = 0.5
        self.far_cut_depth = 0.5
        self.resolution = np.array([1024,1024],'i')
        self.vars = dict(gas = [],star = [], dm = [])
        self.weight = ['density','star/density']
        self.data_maps = []
        # These are details for the case of a HEALPix map
        self.nside = 32
        self.npix = hp.nside2npix(self.nside)

    def save_FITS(self,name,unit_system=basic_conv):
        if self.pov.split('_')[0] != 'mollweide':
            self.save_FITS_image(name,unit_system=basic_conv)
        else:
            self.save_FITS_healpix(name,unit_system=basic_conv)

    def save_FITS_image(self,name,unit_system=basic_conv):
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
                    # Some fields have specific numerical flags at the end
                    # which do not interfere with the units. If that is 
                    # the case, get rid of that last bit
                    try:
                        numflag = int(field.split('/')[1].split('_')[-1])
                        numflag = True
                    except:
                        numflag = False
                    if numflag:
                        sfrstr = field.split('/')[1].split('_')[0] +'_'+ field.split('/')[1].split('_')[1]
                        code_units = get_code_units(sfrstr)
                    else:
                        code_units = get_code_units(field.split('/')[1])
                    temp_map = YTArray(imap[i],code_units,
                                        registry=self.group.obj.unit_registry)
                    first_unit = True
                    for u in code_units.split('*'):
                        if first_unit:
                            units = unit_system[u]
                            first_unit = False
                            units_check = u
                        else:
                            units += '*'+unit_system[u]
                            units_check +='_'+u
                    if units_check in unit_system:
                        units = unit_system[units_check]
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
        centre = [0.0,0.0]
        cdelt = [dx,dy]
        w.wcs.crpix = 0.5*(np.array(self.resolution)+1)
        w.wcs.crval = centre
        w.wcs.cdelt = cdelt
        w.wcs.ctype = ["linear"]*2
        w.wcs.cunit = [unit_system['code_length']]*2
        self.set_wcs(w)

        #TODO: Add further information such that the FITS file is a standalone.

        # Now, just save to disk
        self.hdulist.writeto(name)
    
    def save_FITS_healpix(self,name,unit_system=basic_conv):
        if os.path.exists(name):
            print('WARNING: Overwritting previous HEALPix projection.')
            os.remove(name)
        from yt import YTArray, YTQuantity
        cols = []
        first = True
        counter = 0
        for datatype, varlist in self.vars.items():
            fields = []
            if len(varlist)>0:
                imap = self.data_maps[counter]
                for f in varlist:
                    fields.append(datatype+'/'+f)
                for i,field in enumerate(fields):
                    # Some fields have specific numerical flags at the end
                    # which do not interfere with the units. If that is 
                    # the case, get rid of that last bit
                    try:
                        numflag = int(field.split('/')[1].split('_')[-1])
                        numflag = True
                    except:
                        numflag = False
                    if numflag:
                        sfrstr = field.split('/')[1].split('_')[0] +'_'+ field.split('/')[1].split('_')[1]
                        code_units = get_code_units(sfrstr)
                    else:
                        code_units = get_code_units(field.split('/')[1])
                    temp_map = YTArray(imap[i][0],code_units,
                                        registry=self.group.obj.unit_registry)
                    first_unit = True
                    for u in code_units.split('*'):
                        if first_unit:
                            units = unit_system[u]
                            first_unit = False
                            units_check = u
                        else:
                            units += '*'+unit_system[u]
                            units_check +='_'+u
                    if units_check in unit_system:
                        units = unit_system[units_check]
                    mm = np.array(temp_map.in_units(units))
                    ft = self.getformat(mm[0])
                    cols.append(
                        fits.Column(name=field, format="%s" % ft, 
                                    array=mm, unit=re.sub('()', '', units))
                    )
            counter += 1
        tbhdu = fits.BinTableHDU.from_columns(cols)
        # Add needed keywords
        tbhdu.header["PIXTYPE"] = ("HEALPIX", "HEALPIX pixelisation")
        tbhdu.header["ORDERING"] = ("RING","Pixel ordering scheme, either RING or NESTED")
        tbhdu.header["EXTNAME"] = ("xtension", "Name of this binary table extension")
        tbhdu.header["NSIDE"] = (self.nside, "Resolution parameter of HEALPIX")
        tbhdu.header["COORDSYS"] = ('C',"Ecliptic, Galactic or Celestial (equatorial)")
        tbhdu.header["REDSHIFT"] = (self.group.obj.simulation.redshift, "Snapshot redshift")

        #TODO: Add further information such that the FITS file is a standalone.

        # Now, just save to disk
        tbhdu.writeto(name)

    def set_wcs(self,wcs):
        """Set the WCS coordinate information for all images in a Projection class."""
        wcs.wcs.name = self.pov
        h = wcs.to_header()
        self.wcs = wcs
        for img in self.hdulist:
            for k,v in h.items():
                kk = k
                img.header[kk] = v
    
    def getformat(self,t):
        """Get the FITS convention format string of data type t.
        Parameters
        ----------
        t : data type
        The data type for which the FITS type is requested
        Returns
        -------
        fits_type : str or None
        The FITS string code describing the data type, or None if unknown type.

        TAKEN FROM HEALPY/fitsfuncs.py
        https://github.com/healpy/healpy/blob/d6da2311482de429e08cd031102588d6ade401c8/healpy/fitsfunc.py#L708
        """
        conv = {
            np.dtype(np.bool_): "L",
            np.dtype(np.uint8): "B",
            np.dtype(np.int16): "I",
            np.dtype(np.int32): "J",
            np.dtype(np.int64): "K",
            np.dtype(np.float32): "E",
            np.dtype(np.float64): "D",
            np.dtype(np.complex64): "C",
            np.dtype(np.complex128): "M",
        }
        try:
            if t in conv:
                return conv[t]
        except:
            pass
        try:
            if np.dtype(t) in conv:
                return conv[np.dtype(t)]
        except:
            pass
        try:
            if np.dtype(type(t)) in conv:
                return conv[np.dtype(type(t))]
        except:
            pass
        try:
            if np.dtype(type(t[0])) in conv:
                return conv[np.dtype(type(t[0]))]
        except:
            pass
        try:
            if t is str:
                return "A"
        except:
            pass
        try:
            if type(t) is str:
                return "A%d" % (len(t))
        except:
            pass
        try:
            if type(t[0]) is str:
                l = max(len(s) for s in t)
                return "A%d" % (l)
        except:
            pass

def do_projection(group,vars,weight=['gas/density','star/cumulative'],map_max_size=1024,pov='faceon',lmax=0,lmin=1,window=0.0):
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
            if var_name.split('_')[0] == 'sfr':
                if len(var_name.split('_')) == 3:
                    sfr_name = var_name.split('_')[0] +'_'+var_name.split('_')[1]
                else:
                    sfr_name = var_name.split('_')[0]
                if sfr_name in particle_variables:
                    proj.vars['star'].append(var_name)
                else:
                    raise KeyError('This star variable is not supported. Please check!')
            else:
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
        if weight_type == 'gas':
            if weight_name in common_variables or weight_name in grid_variables:
                proj.weight[0] = weight_name
            else:
                raise KeyError('This gas variable is not supported. Please check!: %s', var)
        elif weight_type == 'star':
            if weight_name in common_variables or weight_name in particle_variables:
                proj.weight[1] = w
            else:
                raise KeyError('This star variable is not supported. Please check!')
        elif weight_type == 'dm':
            if weight_name in common_variables or weight_name in particle_variables:
                proj.weight[1] = w
            else:
                raise KeyError('This DM variable is not supported. Please check!')

    # Setup camera details for the requested POV (Point of View)
    if window != 0.0:
        window = YTArray(window,'kpc',registry=group.obj.unit_registry).in_units('code_length')
    else:
        window = 0.2*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
    centre = vectors.vector()
    centre.x, centre.y, centre.z = group.position[0], group.position[1], group.position[2]
    bulk = vectors.vector()
    velocity = YTArray(group.velocity,'km/s',registry=group.obj.unit_registry).in_units('code_velocity')
    bulk.x, bulk.y, bulk.z = velocity.d[0], velocity.d[1], velocity.d[2]
    if proj.pov == 'faceon':
        axis = vectors.vector()
        norm_L = group.angular_mom['total'].d/np.linalg.norm(group.angular_mom['total'].d)
        up = cartesian_basis['x'] - np.dot(cartesian_basis['x'],norm_L)*norm_L
        up /= np.linalg.norm(up)
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        up_vector = vectors.vector()
        up_vector.x, up_vector.y, up_vector.z = up[0], up[1], up[2]
        rmax = window
        region_size = np.array([2.0*rmax,2.0*rmax],order='F',dtype=np.float64)
        distance = rmax
        far_cut_depth = rmax
    elif proj.pov == 'edgeon':
        axis = vectors.vector()
        norm_L = group.angular_mom['total'].d/np.linalg.norm(group.angular_mom['total'].d)
        los = cartesian_basis['x'] - np.dot(cartesian_basis['x'],norm_L)*norm_L
        los /= np.linalg.norm(los)
        axis.x,axis.y,axis.z = los[0], los[1], los[2]
        up_vector = vectors.vector()
        up_vector.x,up_vector.y,up_vector.z = norm_L[0], norm_L[1], norm_L[2]
        rmax = window
        region_size = np.array([2.0*rmax,2.0*rmax],order='F',dtype=np.float64)
        distance = rmax
        far_cut_depth = rmax
    elif proj.pov == 'x':
        axis = vectors.vector()
        axis.x,axis.y,axis.z = 1.0, 0.0, 0.0
        up_vector = vectors.vector()
        up_vector.x, up_vector.y, up_vector.z = 0.0, 1.0, 0.0
        rmax = window
        region_size = np.array([2.0*rmax,2.0*rmax],order='F',dtype=np.float64)
        distance = rmax
        far_cut_depth = rmax
    elif proj.pov == 'y':
        axis = vectors.vector()
        axis.x,axis.y,axis.z = 0.0, 1.0, 0.0
        up_vector = vectors.vector()
        up_vector.x, up_vector.y, up_vector.z = 1.0, 0.0, 0.0
        rmax = window
        region_size = np.array([2.0*rmax,2.0*rmax],order='F',dtype=np.float64)
        distance = rmax
        far_cut_depth = rmax
    elif proj.pov == 'z':
        axis = vectors.vector()
        axis.x,axis.y,axis.z = 0.0, 0.0, 1.0
        up_vector = vectors.vector()
        up_vector.x, up_vector.y, up_vector.z = 1.0, 0.0, 0.0
        rmax = window
        region_size = np.array([2.0*rmax,2.0*rmax],order='F',dtype=np.float64)
        distance = rmax
        far_cut_depth = rmax
    else:
        print("This point of view is not supported!")
        print("Falling back to 'faceon' (default).")
        axis = vectors.vector()
        norm_L = group.angular_mom['total']/np.linalg.norm(group.angular_mom['total'])
        up = cartesian_basis['y'] - np.dot(cartesian_basis['y'],norm_L)*norm_L
        up /= np.linalg.norm(up)
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        up_vector = vectors.vector()
        up_vector.x, up_vector.y, up_vector.z = up[0], up[1], up[2]
        rmax = window
        region_size = np.array([2.0*rmax,2.0*rmax],order='F')
        distance = rmax
        far_cut_depth = rmax
    
    cam = obs_instruments.init_camera(centre,axis,up_vector,region_size,distance,far_cut_depth,map_max_size-1)

    # Update projection details with the camera ones
    proj.width = region_size
    proj.centre = np.array([group.position[0], group.position[1], group.position[2]])
    proj.distance = distance
    proj.far_cut_depth = far_cut_depth
    proj.los_axis = np.array([cam.los_axis.x,cam.los_axis.y,cam.los_axis.z])
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
    if len(proj.vars['gas']) != 0:
        if lmax != 0:
            maps.projection_hydro(group.obj.simulation.fullpath,cam,bulk,hydro_handler,int(lmax),int(lmin))
        else:
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

    for i in range(0, len(proj.vars['star'])):
        tempstr = 'star/'+proj.vars['star'][i]
        parts_handler.varnames.T.view('S128')[i] = tempstr.ljust(128)
    for i in range(len(proj.vars['star']), len(proj.vars['star'])+len(proj.vars['dm'])):
        tempstr = 'dm/'+proj.vars['dm'][i-len(proj.vars['star'])]
        parts_handler.varnames.T.view('S128')[i] = tempstr.ljust(128)

    # COMPUTE PARTICLES PROJECTION
    if len(proj.vars['star'])+len(proj.vars['dm']) != 0:
        maps.projection_parts(group.obj.simulation.fullpath,cam,bulk,parts_handler)
        # TODO: Weird issue when the direct toto array is given.
        data = np.copy(parts_handler.toto)
        data_star = data[:len(proj.vars['star'])].reshape(len(proj.vars['star']),data.shape[1],data.shape[2])
        proj.data_maps.append(data_star)
        data_dm = data[len(proj.vars['star']):].reshape(len(proj.vars['dm']),data.shape[1],data.shape[2])
        proj.data_maps.append(data_dm)

    return proj

# def plot_single_galaxy_projection(proj_FITS,fields,logscale=True,scalebar=True,redshift=True):
#     """This function uses the projection information in a FITS file following the 
#         OZY format and plots it following the OZY standards."""
    
#     # Make required imports
#     import matplotlib
#     import matplotlib.pyplot as plt
#     import matplotlib.font_manager as fm
#     from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
#     from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
#     from matplotlib.colors import LogNorm
#     import seaborn as sns
#     sns.set(style="white")
#     plt.rc('text', usetex=True)
#     plt.rc('font', family='serif')
#     hfont = {'fontname':'Helvetica'}
#     matplotlib.rc('text', usetex = True)
#     matplotlib.rc('font', **{'family' : "serif"})
#     params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
#     matplotlib.rcParams.update(params)


#     # First,check that FITS file actually exists
#     if not os.path.exists(proj_FITS):
#         raise ImportError('File not found. Please check!')
    
#     # Load FITS file
#     hdul = fits.open(proj_FITS)
#     hdul_fields = [h.header['btype'] for h in hdul]

#     # Check that the required fields for plotting are in this FITS
#     for i,f in enumerate(fields):
#         if f not in hdul_fields:
#             print('The field %s is not included in this file. Ignoring...'%f)
#             del fields[i]
#     if len(fields) == 0:
#         print('Not a single field of the ones you provided are here... Check!')
#         exit

#     # Since everything is fine, we begin plotting…
#     ncolumns = len(fields)
#     fig = plt.figure(figsize=(11,5))
#     grid = AxesGrid(fig, (111),
#                     nrows_ncols = (1, ncolumns),
#                     axes_pad = 0.0,
#                     label_mode = "L",
#                     share_all = True,
#                     cbar_location="top",
#                     cbar_mode="edge",
#                     cbar_size="2%",
#                     cbar_pad=0.0)
#     width_x =  hdul[0].header['CDELT1']*hdul[0].header['NAXIS1']
#     width_y =  hdul[0].header['CDELT2']*hdul[0].header['NAXIS2']
#     ex = [-0.5*width_x,0.5*width_x,-0.5*width_y,0.5*width_y]

#     stellar = False
#     for i in range(0, len(fields)):
#         h = [j for j in range(0,len(hdul)) if hdul[j].header['btype']==fields[i]][0]
#         ax = grid[i].axes
#         ax.axes.xaxis.set_visible(False)
#         ax.axes.yaxis.set_visible(False)
#         if fields[i].split('/')[0] == 'star' or fields[i].split('/')[0] == 'dm':
#             plotting_def = plotting_dictionary[fields[i].split('/')[0]+'_'+fields[i].split('/')[1]]
#             stellar = True
#         else:
#             plotting_def = plotting_dictionary[fields[i].split('/')[1]]
#         if logscale:
#             print(fields[i],np.min(hdul[h].data.T),np.max(hdul[h].data.T))
#             plot = ax.imshow(hdul[h].data.T, cmap=plotting_def['cmap'],
#                             origin='upper',norm=LogNorm(vmin=plotting_def['vmin'],
#                             vmax=plotting_def['vmax']),extent=ex,
#                             interpolation='nearest')
#         else:
#             plot = ax.imshow(hdul[h].data.T, cmap=plotting_def['cmap'],
#                             origin='upper',extent=ex,interpolation='nearest',
#                             vmin=plotting_def['vmin'],vmax=plotting_def['vmax'])
#         fontprops = fm.FontProperties(size=14)
#         if scalebar:
#             scalebar = AnchoredSizeBar(ax.transData,
#                                         3, '3 kpc', 'lower left', 
#                                         pad=0.1,
#                                         color=plotting_def['text_over'],
#                                         frameon=False,
#                                         size_vertical=0.2,
#                                         fontproperties=fontprops)
#             ax.add_artist(scalebar)
#         axcb = fig.colorbar(plot, cax = grid.cbar_axes[i], orientation='horizontal')
#         axcb.set_label(plotting_def['label'], fontsize=16,labelpad=-50, y=0.85)
#         axcb.ax.xaxis.set_ticks_position("top")
#         if redshift:
#             ax.text(0.7, 0.12, 'z = '+str(round(hdul[h].header['redshift'], 2)),
#                     transform=ax.transAxes, fontsize=18,verticalalignment='top',
#                     color=plotting_def['text_over'])

#     fig.subplots_adjust(hspace=0, wspace=0,top = 0.95,bottom = 0.02,left = 0.02,right = 0.98)
#     if stellar:
#         fig.savefig(proj_FITS.split('.')[0]+'_stars.png',format='png',dpi=300)
#     else:
#         fig.savefig(proj_FITS.split('.')[0]+'.png',format='png',dpi=300)

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
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import seaborn as sns
    sns.set(style="dark")
    plt.rcParams["axes.axisbelow"] = False
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    # hfont = {'fontname':'Helvetica'}
    # matplotlib.rc('text', usetex = True)
    # matplotlib.rc('font', **{'family' : "serif"})
    # params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
    # matplotlib.rcParams.update(params)


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
            del fields[i]
    if len(fields) == 0:
        print('Not a single field of the ones you provided are here... Check!')
        exit

    # Since everything is fine, we begin plotting…
    ncol = int(len(fields)/2)
    figsize = plt.figaspect(float(5 * 2) / float(5 * ncol))
    fig = plt.figure(figsize=figsize, facecolor='k', edgecolor='k')
    plot_grid = fig.add_gridspec(2, ncol, wspace=0, hspace=0,left=0,right=1, bottom=0, top=1)
    ax = []
    for i in range(0,2):
        ax.append([])
        for j in range(0,ncol):
            ax[i].append(fig.add_subplot(plot_grid[i,j]))
    ax = np.asarray(ax)

    width_x =  hdul[0].header['CDELT1']*hdul[0].header['NAXIS1']
    width_y =  hdul[0].header['CDELT2']*hdul[0].header['NAXIS2']
    ex = [-0.5*width_x,0.5*width_x,-0.5*width_y,0.5*width_y]

    stellar = False
    for i in range(0, ax.shape[0]):
        for j in range(0, ax.shape[1]):
            ivar = i*ax.shape[1] + j
            if ivar >= len(fields):
                # Clear that extra panel
                ax[i,j].get_xaxis().set_visible(False)
                ax[i,j].get_yaxis().set_visible(False)
                break
            ax[i,j].set_xlim([-0.5*width_x,0.5*width_y])
            ax[i,j].set_ylim([-0.5*width_x,0.5*width_y])
            ax[i,j].axes.xaxis.set_visible(False)
            ax[i,j].axes.yaxis.set_visible(False)
            ax[i,j].axis('off')
            h = [k for k in range(0,len(hdul)) if hdul[k].header['btype']==fields[ivar]][0]

            if fields[ivar].split('/')[0] == 'star' or fields[ivar].split('/')[0] == 'dm':
                plotting_def = plotting_dictionary[fields[ivar].split('/')[0]+'_'+fields[ivar].split('/')[1]]
                stellar = True
            else:
                plotting_def = plotting_dictionary[fields[ivar].split('/')[1]]
            if logscale:
                print(fields[ivar],np.min(hdul[h].data.T),np.max(hdul[h].data.T))
                plot = ax[i,j].imshow(np.log10(hdul[h].data.T), cmap=plotting_def['cmap'],
                                origin='upper',vmin=np.log10(plotting_def['vmin']),
                                vmax=np.log10(plotting_def['vmax']),extent=ex,
                                interpolation='nearest')
            else:
                plot = ax[i,j].imshow(hdul[h].data.T, cmap=plotting_def['cmap'],
                                origin='upper',extent=ex,interpolation='nearest',
                                vmin=plotting_def['vmin'],vmax=plotting_def['vmax'])

            cbaxes = inset_axes(ax[i,j], width="80%", height="5%", loc='lower center')
            cbar = fig.colorbar(plot, cax=cbaxes, orientation='horizontal')
            if logscale:
                cbar.set_label(plotting_def['label_log'],color=plotting_def['text_over'],fontsize=20,labelpad=-25, y=0.85,weight='bold')
            else:
                cbar.set_label(plotting_def['label'],color=plotting_def['text_over'],fontsize=16,labelpad=-25, y=0.85)
            cbar.ax.xaxis.label.set_font_properties(matplotlib.font_manager.FontProperties(weight='bold',size=8))
            cbar.ax.tick_params(axis='x', pad=-7, labelsize=8,labelcolor=plotting_def['text_over'])
            cbar.ax.tick_params(length=0,width=0)

            if redshift and i==0 and j==0:
                ax[i,j].text(0.05, 0.90, r'$z = ${z:.2f}'.format(z=hdul[h].header['redshift']), # Redshift
                                    verticalalignment='bottom', horizontalalignment='left',
                                    transform=ax[i,j].transAxes,
                                    color=plotting_def['text_over'], fontsize=10,fontweight='bold')

            fontprops = fm.FontProperties(size=10,weight='bold')
            if scalebar and i==0 and j==0:
                scalebar = AnchoredSizeBar(ax[i,j].transData,
                                            1, '1 kpc', 'upper right', 
                                            pad=0.1,
                                            color=plotting_def['text_over'],
                                            frameon=False,
                                            size_vertical=0.1,
                                            fontproperties=fontprops)
                ax[i,j].add_artist(scalebar)


    fig.subplots_adjust(hspace=0,wspace=0,left=0,right=1, bottom=0, top=1)
    if stellar:
        fig.savefig(proj_FITS.split('.')[0]+'_stars.png',format='png',dpi=300)
    else:
        fig.savefig(proj_FITS.split('.')[0]+'.png',format='png',dpi=300)

def plot_lupton_rgb_projection(proj_FITS,fields,stars=False,scalebar=True,redshift=True):
    """This function uses the projection information in a FITS file following the 
        OZY format and combines three variables into an RGB image using the Lupton et al.
        (2004) method."""
    
    # Make required imports
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from matplotlib.colors import LogNorm
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from astropy.visualization import make_lupton_rgb
    import matplotlib.patches as mpatches
    from matplotlib.colors import Normalize
    import seaborn as sns
    sns.set(style="dark")
    plt.rcParams["axes.axisbelow"] = False

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
            del fields[i]
    if len(fields) == 0:
        print('Not a single field of the ones you provided are here... Check!')
        exit
    if len(fields) != 3:
        print('You need to provide three different fields to make RGB images!')
        exit
    
    # Since everything is fine, we begin plotting…
    figsize = plt.figaspect(float(10) / float(10))
    fig, ax = plt.subplots(1,1,figsize=figsize, facecolor='k', edgecolor='k')

    width_x =  hdul[0].header['CDELT1']*hdul[0].header['NAXIS1']
    width_y =  hdul[0].header['CDELT2']*hdul[0].header['NAXIS2']
    ex = [-0.5*width_x,0.5*width_x,-0.5*width_y,0.5*width_y]

    rgb_colours = ['red','green','blue']
    images = []
    custom_lines = []

    for i in range(0, 3):
        ax.set_xlim([-0.5*width_x,0.5*width_y])
        ax.set_ylim([-0.5*width_x,0.5*width_y])
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
        ax.axis('off')
        h = [k for k in range(0,len(hdul)) if hdul[k].header['btype']==fields[i]][0]

        if fields[i].split('/')[0] == 'star' or fields[i].split('/')[0] == 'dm':
            plotting_def = plotting_dictionary[fields[i].split('/')[0]+'_'+fields[i].split('/')[1]]
        else:
            plotting_def = plotting_dictionary[fields[i].split('/')[1]]

        data = hdul[h].data.T
        print(fields[i],data.min(),data.max())
        data[data < plotting_def['vmin_galaxy']] = plotting_def['vmin_galaxy']
        data[data > plotting_def['vmax_galaxy']] = plotting_def['vmax_galaxy']
        data = (np.log10(data)-np.log10(plotting_def['vmin_galaxy']*0.99))/(np.log10(plotting_def['vmax_galaxy']) - np.log10(plotting_def['vmin_galaxy']*0.99))

        images.append(data)
        custom_lines.append(mpatches.Patch(color=rgb_colours[i],label=plotting_def['label']))
    
    rgb_default = make_lupton_rgb(images[0], images[1], images[2],Q=10, stretch=0.5)

    ax.imshow(rgb_default,origin='upper',extent=ex,interpolation='nearest')

    if stars:
        h = [k for k in range(0,len(hdul)) if hdul[k].header['btype']=='star/mass'][0]
        stars = np.log10(hdul[h].data.T)
        plotting_def = plotting_dictionary['star_mass']
        ax.imshow(stars,origin='upper',cmap=plotting_def['cmap'],
                    extent=ex,interpolation='nearest',
                    vmin=np.log10(plotting_def['vmin']),
                    vmax=np.log10(plotting_def['vmax']),alpha=0.4)

    if redshift:
        ax.text(0.05, 0.90, r'$z = ${z:.2f}'.format(z=hdul[h].header['redshift']), # Redshift
                    verticalalignment='bottom', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='white', fontsize=16,fontweight='bold')

    fontprops = fm.FontProperties(size=16,weight='bold')
    if scalebar:
        scalebar = AnchoredSizeBar(ax.transData,
                                    1, '1 kpc', 'upper right', 
                                    pad=0.1,
                                    color='white',
                                    frameon=False,
                                    size_vertical=0.1,
                                    fontproperties=fontprops)
        ax.add_artist(scalebar)
    ax.legend(handles=custom_lines,loc='lower right',frameon=False,fontsize=12,labelcolor='white')
    fig.subplots_adjust(hspace=0,wspace=0,left=0,right=1, bottom=0, top=1)
    fig.savefig(proj_FITS.split('.')[0]+'_rgb.png',format='png',dpi=300)


def do_healpix_projection(group,vars,weight=['gas/density','star/age'],nside=32,pov='edgeon'):
    """Function which computes a 2D spherical projection of particular object using the HEALPix
        pixelisation scheme."""
    
    if not isinstance(weight,list):
        weight = [weight]
    if len(weight)>2:
        print('Only just one weight per data type (i.e. grid and particles)!')
        exit
    proj = Projection(group)
    proj.pov = 'mollweide_'+str(pov)

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
        if weight_type == 'gas':
            if weight_name in common_variables or weight_name in grid_variables:
                proj.weight[0] = weight_name
            else:
                raise KeyError('This gas variable is not supported. Please check!: %s', var)
        elif weight_type == 'star':
            if weight_name in common_variables or weight_name in particle_variables:
                proj.weight[1] = w
            else:
                raise KeyError('This star variable is not supported. Please check!')
        elif weight_type == 'dm':
            if weight_name in common_variables or weight_name in particle_variables:
                proj.weight[1] = w
            else:
                raise KeyError('This DM variable is not supported. Please check!')

    # Setup region details
    if pov == 'edgeon':
        reg = geometrical_regions.region()
        reg.name = 'sphere'
        centre = vectors.vector()
        centre.x, centre.y, centre.z = group.position[0], group.position[1], group.position[2]
        reg.centre = centre
        bulk = vectors.vector()
        velocity = YTArray(group.velocity,'km/s',registry=group.obj.unit_registry).in_units('code_velocity')
        bulk.x, bulk.y, bulk.z = velocity.d[0], velocity.d[1], velocity.d[2]
        reg.bulk_velocity = bulk
        norm_L = group.angular_mom['total'].d/np.linalg.norm(group.angular_mom['total'].d)
        axis = vectors.vector()
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        reg.axis = axis
        rmax = 0.4*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        reg.rmin = 0.85*rmax
        reg.rmax = 1.15*rmax

        # Update projection details with the ones used for the region
        proj.up_vector = norm_L
        proj.centre = group.position[:]
        proj.nside = nside
    else:
        print('This POV for a HEALPix projection is not supported...')
        print('Please check!')
        exit
    
    # Create projection_handler Fortran derived type for the results of the hydro data projection
    hydro_handler = maps.projection_handler()
    hydro_handler.type = pov
    hydro_handler.nvars = len(proj.vars['gas'])
    hydro_handler.weightvar = proj.weight[0]
    maps.allocate_projection_handler(hydro_handler)
    for i in range(0, len(proj.vars['gas'])):
        hydro_handler.varnames.T.view('S128')[i] = proj.vars['gas'][i].ljust(128)
    
    # COMPUTE HYDRO PROJECTION
    maps.healpix_hydro(group.obj.simulation.fullpath,reg,nside,hydro_handler)
    # TODO: Weird issue when the direct toto array is given.
    data = np.copy(hydro_handler.toto)
    proj.data_maps.append(data)

    # TODO: Add particle projections

    return proj
    
def plot_single_galaxy_healpix(proj_FITS,fields,logscale=True,redshift=False):
    """This function uses the HEALPix projection information in a FITS file following
        the OZY format and plots it using the OZY standards."""

    # Make required imports
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from matplotlib.colors import LogNorm,SymLogNorm
    import seaborn as sns
    from reproject import reproject_from_healpix
    from astropy.wcs import WCS
    from astropy.visualization.wcsaxes.frame import EllipticalFrame
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    sns.set(style="white")
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    hfont = {'fontname':'Helvetica'}
    matplotlib.rc('text', usetex = True)
    matplotlib.rc('font', **{'family' : "serif"})
    params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
    matplotlib.rcParams.update(params)

    # First, check that the FITS file actually exists
    if not os.path.exists(proj_FITS):
        raise ImportError('File not found. Please check!')
    
    # Load FITS file
    hdul = fits.open(proj_FITS)[1]
    hdul_fields = hdul.data.columns.names

    # Check that the required fields for plotting are in this FITS
    for i,f in enumerate(fields):
        if f not in hdul_fields:
            print('The field %s is not included in this file. Ignoring...'%f)
            del f[i]
    if len(fields) == 0:
        print('Not a single field of the ones you provided are here... Check!')
        exit
    
    # Make adaptations for HEALPix and Astropy
    target_wcs = WCS(target_header)
    print(target_wcs)
    # Since everything is fine, we begin plotting…
    nside = hdul.header["NSIDE"]
    ncolumns = int(len(fields)/2)
    figsize = plt.figaspect(float(7 * 2) / float(11 * ncolumns))
    fig = plt.figure(figsize=figsize)
    plot_grid = fig.add_gridspec(2, ncolumns)
    axes = []
    for i in range(0,2):
        axes.append([])
        for j in range(0,ncolumns):
            axes[i].append(fig.add_subplot(plot_grid[i,j],
                            projection=target_wcs,
                            frame_class=EllipticalFrame))
    axes = np.asarray(axes)

    stellar  = False
    for i in range(0, axes.shape[0]):
        for j in range(0, axes.shape[1]):
            ivar = i*axes.shape[1] + j
            ax = axes[i,j]
            if ivar >= len(fields) :
                # Clear that extra panel
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                break
            data = hdul.data[fields[ivar]]
            invalid_index = np.argwhere(np.isnan(data))
            theta,phi = hp.pix2ang(nside,invalid_index)
            data[invalid_index] = np.nanmin(data)
            array, footprint = reproject_from_healpix((data, 'C'),
                                                        target_header,nested=False)
            print(fields[ivar],np.nanmin(array),np.nanmax(array))
            
            if fields[ivar].split('/')[0] == 'star' or fields[ivar].split('/')[0] == 'dm':
                plotting_def = plotting_dictionary[fields[ivar].split('/')[0]+'_'+fields[ivar].split('/')[1]]
                stellar = True
            else:
                plotting_def = plotting_dictionary[fields[ivar].split('/')[1]]
            if logscale:
                if fields[ivar].split('/')[1].split('_')[0] != 'v':
                    if fields[ivar].split('/')[1] == 'density':
                        vmin = plotting_def['vmin']/100
                        vmax = plotting_def['vmax']/100
                    elif fields[ivar].split('/')[1] == 'temperature':
                        vmin = plotting_def['vmin']*100
                        vmax = plotting_def['vmax']*10
                    elif fields[ivar].split('/')[1] == 'xHII':
                        vmin = plotting_def['vmin']*100
                        vmax = plotting_def['vmax']
                    # elif fields[ivar].split('/')[1] == 'magnetic_energy_density':
                    #     vmin = plotting_def['vmin']/1e12
                    #     vmax = plotting_def['vmax']*10
                    else:
                        vmin = plotting_def['vmin']
                        vmax = plotting_def['vmax']
                    plot = ax.imshow(array, cmap=plotting_def['cmap'],
                                    norm=LogNorm(vmin=vmin,
                                    vmax=vmax),
                                    interpolation='nearest')
                else:
                    plot = ax.imshow(array, cmap=plotting_def['cmap'],
                                    vmin=plotting_def['vmin'],
                                    vmax=plotting_def['vmax'],
                                    interpolation='nearest')
            else:
                plot = ax.imshow(array, cmap=plotting_def['cmap'],
                                interpolation='nearest',
                                vmin=plotting_def['vmin'],vmax=plotting_def['vmax'])
            ax.coords.grid(color=plotting_def['text_over'],linewidth=0.8,alpha=0.4,linestyle=':')
            # ax.coords['ra'].set_ticklabel(color=plotting_def['text_over'])
            ax.coords['ra'].set_ticklabel_visible(False)
            ax.coords['ra'].set_ticks_visible(False)
            ax.coords['dec'].set_ticks_position('lb')
            ax.coords['dec'].set_ticklabel_position('lb')
            fontprops = fm.FontProperties(size=14)
            if i == 0:
                cbaxes = inset_axes(ax, width="80%", height="7%", loc='upper center',
                                    bbox_to_anchor=(0.0, 0.2, 1, 1.),
                                    bbox_transform=ax.transAxes)
                axcb = fig.colorbar(plot, cax = cbaxes, orientation='horizontal',pad=2)
                axcb.set_label(plotting_def['label'], fontsize=16,labelpad=-50, y=0.85)
                axcb.ax.xaxis.set_ticks_position("top")
            else:
                cbaxes = inset_axes(ax, width="80%", height="7%", loc='lower center',
                                    bbox_to_anchor=(0.0, -0.2, 1., 1.),
                                    bbox_transform=ax.transAxes)
                axcb = fig.colorbar(plot, cax = cbaxes, orientation='horizontal',pad=0.3)
                axcb.set_label(plotting_def['label'], fontsize=16)
    if redshift:
        fig.text(0.5, 0.5, 'z = '+str(round(hdul.header['REDSHIFT'], 2)),
                fontsize=18,va='center',
                color=plotting_def['text_over'])
    # fig.subplots_adjust(hspace=0, wspace=0,top = 0.95,bottom = 0.02,left = 0.02,right = 0.98)
    if stellar:
        fig.savefig(proj_FITS.split('.')[0]+'_stars.png',format='png',dpi=300)
    else:
        fig.savefig(proj_FITS.split('.')[0]+'.png',format='png',dpi=300)