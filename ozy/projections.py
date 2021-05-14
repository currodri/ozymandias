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
from projections import obs_instruments
from projections import maps
from projections import vectors
from projections import geometrical_regions
from astropy.io import fits
from astropy.wcs import WCS


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
                    code_units = get_code_units(field.split('/')[1])
                    print(fields[i],code_units,np.min(imap[i]),np.max(imap[i]))
                    temp_map = YTArray(imap[i],code_units,
                                        registry=self.group.obj.unit_registry)
                    first_unit = True
                    print(code_units.split('*'))
                    for u in code_units.split('*'):
                        print(u)
                        if first_unit:
                            units = unit_system[u]
                            first_unit = False
                            units_check = u
                        else:
                            units += '*'+unit_system[u]
                            units_check +='_'+u
                    print(units_check)
                    if units_check in unit_system:
                        units = unit_system[units_check]
                    print(units,np.array(temp_map.in_units(units)).min(),np.array(temp_map.in_units(units)).max())
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
                    code_units = get_code_units(field.split('/')[1])
                    temp_map = YTArray(imap[i][0],code_units,
                                        registry=self.group.obj.unit_registry)
                    first_unit = True
                    for u in code_units.split('*'):
                        if first_unit:
                            units = unit_system[u]
                            first_unit = False
                        else:
                            units += '*'+unit_system[u]
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
        up = cartesian_basis['y'] - np.dot(cartesian_basis['y'],norm_L)*norm_L
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
            del fields[i]
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
            print(fields[i],np.min(hdul[h].data.T),np.max(hdul[h].data.T))
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
        rmax = 0.2*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        reg.rmin = 0.9*rmax
        reg.rmax = 1.1*rmax

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
            array, footprint = reproject_from_healpix((data, 'C'),
                                                        target_header,nested=False)
            if fields[ivar].split('/')[0] == 'star' or fields[ivar].split('/')[0] == 'dm':
                plotting_def = plotting_dictionary[fields[ivar].split('/')[0]+'_'+fields[ivar].split('/')[1]]
                stellar = True
            else:
                plotting_def = plotting_dictionary[fields[ivar].split('/')[1]]
            if logscale:
                if fields[ivar].split('/')[1].split('_')[0] != 'v':
                    if fields[ivar].split('/')[1] == 'density':
                        vmin = plotting_def['vmin']/10
                        vmax = plotting_def['vmax']/100
                    elif fields[ivar].split('/')[1] == 'temperature':
                        vmin = plotting_def['vmin']*100
                        vmax = plotting_def['vmax']*10
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