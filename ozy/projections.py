import numpy as np
import h5py
import os
import ozy
from ozy.utils import init_filter, init_region, structure_regions, check_need_neighbours, get_code_units, get_plotting_def
from variables_settings import geometrical_variables,raw_gas_variables,\
    raw_star_variables,raw_dm_variables,derived_gas_variables,\
    derived_star_variables,derived_dm_variables,gravity_variables,\
    basic_conv,circle_dictionary
import re
from unyt import unyt_quantity,unyt_array
import healpy as hp
from astropy.io import fits
from astropy.wcs import WCS
from vis import obs_instruments
from vis import maps
from vis import vectors
from vis import io_ramses

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
        self.obj = group.obj
        self.pov = 'faceon'
        self.los_axis = np.array([0.,0.,1.])
        self.up_vector = np.array([0.,1.,0.])
        self.centre = np.array([0.5,0.5,0.5])
        self.width = np.array([1.,1.])
        self.distance = 0.5
        self.far_cut_depth = 0.5
        self.resolution = np.array([1024,1024],'i')
        self.vars = dict(gas = [],star = [], dm = [])
        self.weight = dict(gas = [],star = [], dm = [])
        self.data_maps = []
        self.filters = []
        # These are details for the case of a HEALPix map
        self.nside = 32
        self.npix = hp.nside2npix(self.nside)

    def save_FITS(self,name,unit_system=basic_conv):
        if self.pov.split('_')[0] != 'mollweide':
            self.save_FITS_image(name,unit_system=unit_system)
        else:
            self.save_FITS_healpix(name,unit_system=unit_system)

    def save_FITS_image(self,name,unit_system=basic_conv):
        if os.path.exists(name):
            print('WARNING: Overwritting previous projection.')
            os.remove(name)
        self.hdulist = fits.HDUList()
        first = True
        nfilter = len(self.filters)
        for nf in range(0,nfilter):
            filter_name = self.filters[nf]['name']
            counter = 0
            for datatype, varlist in self.vars.items():
                fields = []
                if len(varlist)>0:
                    imap = self.data_maps[counter][nf]
                    if not isinstance(imap,list) and not isinstance(imap,np.ndarray):
                        imap = [imap]
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
                            plotting_def = get_plotting_def(sfrstr)
                            units = plotting_def['units']
                        else:
                            code_units = get_code_units(field.split('/')[1])
                            if datatype == 'star' or datatype == 'dm':
                                field_str = field.split('/')[0] + '_' + field.split('/')[1]
                            else:
                                field_str = field.split('/')[1]
                            plotting_def = get_plotting_def(field_str)
                            units = plotting_def['units']                            
                        temp_map = self.group.obj.array(imap[i],code_units)
                        # first_unit = True
                        # make_div = False
                        # if len(code_units.split('/')) != 1:
                        #     make_div = True
                        #     code_units = code_units.replace('/','*')
                        # for u in code_units.split('*'):
                        #     if first_unit:
                        #         units = unit_system[u]
                        #         first_unit = False
                        #         units_check = u
                        #     else:
                        #         if make_div:
                        #             units += '/'+unit_system[u]
                        #         else:
                        #             units += '*'+unit_system[u]
                        #         units_check +='_'+u
                        # if units_check in unit_system:
                        #     units = unit_system[units_check]
                        
                        if first:
                            hdu = fits.PrimaryHDU(np.array(temp_map.in_units(units)))
                            first = False
                        else:
                            hdu = fits.ImageHDU(np.array(temp_map.in_units(units)))
                        if filter_name != 'none':
                            hdu.name = field+'/'+filter_name
                            hdu.header["btype"] = field+'/'+filter_name
                        else:
                            hdu.name = field
                            hdu.header["btype"] = field
                        hdu.header["bunit"] = re.sub('()', '', units)
                        hdu.header["redshift"] = self.group.obj.simulation.redshift
                        hdu.header["los_x"] = self.los_axis[0]
                        hdu.header["los_y"] = self.los_axis[1]
                        hdu.header["los_z"] = self.los_axis[2]
                        hdu.header["up_x"] = self.up_vector[0]
                        hdu.header["up_y"] = self.up_vector[1]
                        hdu.header["up_z"] = self.up_vector[2]
                        hdu.header["centre_x"] = float(self.centre[0].in_units(unit_system['code_length']).d)
                        hdu.header["centre_y"] = float(self.centre[1].in_units(unit_system['code_length']).d)
                        hdu.header["centre_z"] = float(self.centre[2].in_units(unit_system['code_length']).d)
                        self._save_filterinfo(hdu,self.filters[nf])
                        self.hdulist.append(hdu)
                counter += 1
        # Setup WCS in the coordinate systems of the camera
        w = WCS(header=self.hdulist[0].header,naxis=2)
        dx = self.group.obj.array(self.width[0],'code_length').in_units(unit_system['code_length']).d
        dy = self.group.obj.array(self.width[1],'code_length').in_units(unit_system['code_length']).d
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
        cols = []
        first = True
        nfilter = len(self.filters)
        for nf in range(0, nfilter):
            filter_name = self.filters[nf]['name']
            counter = 0
            for datatype, varlist in self.vars.items():
                fields = []
                if len(varlist)>0:
                    imap = self.data_maps[counter][nf]
                    if not isinstance(imap,list) and not isinstance(imap,np.ndarray):
                        imap = [imap]
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
                        temp_map = self.group.obj.array(imap[i][0],code_units)
                        make_div = False
                        first_unit = True
                        if len(code_units.split('/')) != 1:
                            make_div = True
                            code_units = code_units.replace('/','*')
                        for u in code_units.split('*'):
                            if first_unit:
                                units = unit_system[u]
                                first_unit = False
                                units_check = u
                            else:
                                if make_div:
                                    units += '/'+unit_system[u]
                                else:
                                    units += '*'+unit_system[u]
                                units_check +='_'+u
                        if units_check in unit_system:
                            units = unit_system[units_check]
                        mm = np.array(temp_map.in_units(units))
                        ft = self.getformat(mm[0])
                        if filter_name != 'none':
                            cols.append(
                                fits.Column(name=field+'/'+filter_name, format="%s" % ft, 
                                            array=mm, unit=re.sub('()', '', units))
                            )
                        else:
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
        tbhdu.header["los_x"] = self.los_axis[0]
        tbhdu.header["los_y"] = self.los_axis[1]
        tbhdu.header["los_z"] = self.los_axis[2]
        tbhdu.header["up_x"] = self.up_vector[0]
        tbhdu.header["up_y"] = self.up_vector[1]
        tbhdu.header["up_z"] = self.up_vector[2]
        tbhdu.header["centre_x"] = float(self.centre[0].in_units(unit_system['code_length']).d)
        tbhdu.header["centre_y"] = float(self.centre[1].in_units(unit_system['code_length']).d)
        tbhdu.header["centre_z"] = float(self.centre[2].in_units(unit_system['code_length']).d)
        self._save_filterinfo(tbhdu,self.filters[nf])

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

    def _get_python_filter(self,filt):
        """Save the Fortran derived type as a dictionary inside the PhaseDiagram class (only the necessary info)."""
        self.filters.append(dict())
        self.filters[-1]['name'] = filt.name.decode().split(' ')[0]
        self.filters[-1]['conditions'] = []
        if self.filters[-1]['name'] != 'none':
            for i in range(0, filt.ncond):
                particle = False
                cond_var = filt.cond_vars.T.view('S128')[i][0].decode().split(' ')[0]
                if cond_var.split('/')[0] == 'star' or cond_var.split('/')[0] == 'dm':
                    particle = True
                    corrected_part_var = cond_var.split('/')[0] + '_' + cond_var.split('/')[1]
                cond_op = filt.cond_ops.T.view('S2')[i][0].decode().split(' ')[0]
                if particle:
                    cond_units = get_code_units(cond_var.split('/')[1])
                else:
                    cond_units = get_code_units(cond_var)
                cond_value = self.obj.quantity(filt.cond_vals[i], str(cond_units))
                if particle:
                    cond_str = corrected_part_var+'/'+cond_op+'/'+str(cond_value.d)+'/'+cond_units
                else:
                    cond_str = cond_var+'/'+cond_op+'/'+str(cond_value.d)+'/'+cond_units
                self.filters[-1]['conditions'].append(cond_str)
    def _save_filterinfo(self,hdu,filt):
        """"This function unravels the information contained in a filter object into the FITS header"""
        if filt['name'] != 'none':
            for i in range(0, len(filt['conditions'])):
                hdu.header["cond_"+str(i)] = filt['conditions'][i]

def do_projection(group,vars,weight=['gas/density','star/cumulative'],map_max_size=1024,
                    pov='faceon',lmax=100,lmin=1,type_projection='gauss_deposition', window=(0.0,'kpc'),
                    rmin=(0.0,'rvir'), rmax=(0.2,'rvir'), xmin=(0.0,'rvir'), xmax=(0.2,'rvir'),
                    ymin=(0.0,'rvir'), ymax=(0.2,'rvir'),zmin=(0.0,'rvir'), zmax=(0.2,'rvir'),
                    mycentre=([0.5,0.5,0.5],'rvir'), myaxis=np.array([1.,0.,0.]),up_axis=np.array([1.,0.,0.]),
                    thickness=(0.0,'kpc'),
                    nexp_factor = 1.0,
                    tag_file=None,
                    inverse_tag=False, remove_subs=False,
                    filter_conds=['none'],filter_name=['none'],
                    verbose=False):
    """Function which computes a 2D projection centered on an objected from an OZY file."""

    if isinstance(group,ozy.Snapshot):
        obj = group
        group = ozy.group.Group(obj)
        use_snapshot = True
    else:
        obj = group.obj
        use_snapshot = False

    if not isinstance(weight,list):
        weight = [weight]
    if len(weight)>2:
        print('Only just one weight per data type (i.e. grid and particles)!')
        exit
    proj = Projection(group)
    proj.pov = pov
    use_neigh = False
    
    for var in vars:
        var_type = var.split('/')[0]
        var_name = var.split('/')[1]
        if var_type == 'gas':
            use_neigh = check_need_neighbours(var_name) or use_neigh
            if var_name in geometrical_variables or var_name in raw_gas_variables \
                or var_name in derived_gas_variables or var_name in gravity_variables:
                proj.vars['gas'].append(var_name)
            else:
                raise KeyError('This gas variable is not supported. Please check!: %s', var)
        elif var_type == 'star':
            if var_name.split('_')[0] == 'sfr':
                if len(var_name.split('_')) == 3:
                    sfr_name = var_name.split('_')[0] +'_'+var_name.split('_')[1]
                else:
                    sfr_name = var_name.split('_')[0]
                if sfr_name in derived_star_variables:
                    proj.vars['star'].append(var_name)
                else:
                    raise KeyError('This star variable is not supported. Please check!')
            else:
                if var_name in geometrical_variables or var_name in raw_part_variables \
                    or var_name in derived_part_variables or var_name in star_variables:
                    proj.vars['star'].append(var_name)
                else:
                    raise KeyError('This star variable is not supported. Please check!')
        elif var_type == 'dm':
            if var_name in geometrical_variables or var_name in raw_dm_variables \
                or var_name in derived_dm_variables:
                proj.vars['dm'].append(var_name)
            else:
                raise KeyError('This DM variable is not supported. Please check!')

    for w in weight:
        weight_type = w.split('/')[0]
        weight_name = w.split('/')[1]
        if weight_type == 'gas':
            use_neigh = check_need_neighbours(weight_name) or use_neigh
            if weight_name in geometrical_variables or weight_name in raw_gas_variables \
                or weight_name in derived_gas_variables or weight_name in gravity_variables:
                proj.weight['gas'].append(weight_name)
            else:
                raise KeyError('This gas variable is not supported. Please check!: %s', var)
        elif weight_type == 'star':
            if weight_name in geometrical_variables or weight_name in raw_star_variables \
                or weight_name in derived_star_variables:
                proj.weight['star'].append(weight_name)
            else:
                raise KeyError('This star variable is not supported. Please check!')
        elif weight_type == 'dm':
            if weight_name in geometrical_variables or weight_name in raw_dm_variables \
                or weight_name in derived_dm_variables:
                proj.weight['dm'].append(weight_name)
            else:
                raise KeyError('This DM variable is not supported. Please check!')
    if use_neigh and verbose:
        print('At least one variable needs neighbours!')

    # Setup camera details for the requested POV (Point of View)
    boxlen = obj.simulation.boxsize.to('code_length').d   
    if use_snapshot:
        group.position = obj.array(mycentre[0],mycentre[1])
        window = obj.quantity(window[0],window[1]).in_units('code_length').d
        region_axis = vectors.vector()
        region_axis.x,region_axis.y,region_axis.z = 0.,0.,1.
        bulk = vectors.vector()
        
    else:
        if window[0] == 0.0 or window[1] == 'rvir':
            if group.type == 'halo':
                rvir = group.virial_quantities['radius'].d
            else:
                rvir = group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        if window[0] == 0.0:
            if group.type == 'halo':
                window = 1.2*rvir
            else:
                window = 0.2*rvir
        else:
            if window[1] == 'rvir':
                window = window[0]*rvir
            else:
                window = obj.quantity(window[0],window[1]).in_units('code_length').d
        centre = vectors.vector()
        centre.x, centre.y, centre.z = group.position[0], group.position[1], group.position[2]
        bulk = vectors.vector()
        region_axis = vectors.vector()
        norm_L = group.angular_mom['total'].d/np.linalg.norm(group.angular_mom['total'].d)
        region_axis.x,region_axis.y,region_axis.z = norm_L[0], norm_L[1], norm_L[2]

        if group.type != 'halo':
            velocity = group.velocity.in_units('code_velocity')
            bulk.x, bulk.y, bulk.z = velocity.d[0], velocity.d[1], velocity.d[2]

    if proj.pov == 'faceon':
        centre = vectors.vector()
        centre.x, centre.y, centre.z = group.position[0], group.position[1], group.position[2]
        axis = vectors.vector()
        norm_L = group.angular_mom['total'].d/np.linalg.norm(group.angular_mom['total'].d)
        up = cartesian_basis['x'] - np.dot(cartesian_basis['x'],norm_L)*norm_L
        up /= np.linalg.norm(up)
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        up_vector = vectors.vector()
        up_vector.x, up_vector.y, up_vector.z = up[0], up[1], up[2]
        rmax = window
        region_size = np.array([2.0*rmax,2.0*rmax],order='F',dtype=np.float64)
        if thickness[0] != 0.0:
            thickness_value = group.obj.quantity(thickness[0],thickness[1]).in_units('code_length').d
            distance = 0.5*thickness_value
            far_cut_depth = 0.5*thickness_value
        else:
            distance = rmax
            far_cut_depth = rmax
        enclosing_sphere_p = group.position
        enclosing_sphere_r = rmax
    elif proj.pov == 'edgeon':
        centre = vectors.vector()
        centre.x, centre.y, centre.z = group.position[0], group.position[1], group.position[2]
        axis = vectors.vector()
        norm_L = group.angular_mom['total'].d/np.linalg.norm(group.angular_mom['total'].d)
        los = cartesian_basis['x'] - np.dot(cartesian_basis['x'],norm_L)*norm_L
        los /= np.linalg.norm(los)
        axis.x,axis.y,axis.z = -los[0], -los[1], -los[2]
        up_vector = vectors.vector()
        up_vector.x,up_vector.y,up_vector.z = norm_L[0], norm_L[1], norm_L[2]
        rmax = window
        region_size = np.array([2.0*rmax,2.0*rmax],order='F',dtype=np.float64)
        if thickness[0] != 0.0:
            thickness_value = group.obj.quantity(thickness[0],thickness[1]).in_units('code_length').d
            distance = 0.5*thickness_value
            far_cut_depth = 0.5*thickness_value
        else:
            distance = rmax
            far_cut_depth = rmax
        enclosing_sphere_p = group.position
        enclosing_sphere_r = rmax
    elif proj.pov == 'x':
        centre = vectors.vector()
        centre.x, centre.y, centre.z = group.position[0].to('code_length')/boxlen, group.position[1].to('code_length')/boxlen, group.position[2].to('code_length')/boxlen
        axis = vectors.vector()
        axis.x,axis.y,axis.z = 1.0, 0.0, 0.0
        up_vector = vectors.vector()
        up_vector.x, up_vector.y, up_vector.z = 0.0, 1.0, 0.0
        rmax = window
        region_size = np.array([2.0*rmax,2.0*rmax],order='F',dtype=np.float64)
        if thickness[0] != 0.0:
            thickness_value = group.obj.quantity(thickness[0],thickness[1]).in_units('code_length').d
            distance = 0.5*thickness_value
            far_cut_depth = 0.5*thickness_value
        else:
            distance = rmax
            far_cut_depth = rmax
        enclosing_sphere_p = group.position
        enclosing_sphere_r = rmax
    elif proj.pov == 'y':
        centre = vectors.vector()
        centre.x, centre.y, centre.z = group.position[0].to('code_length')/boxlen, group.position[1].to('code_length')/boxlen, group.position[2].to('code_length')/boxlen
        axis = vectors.vector()
        axis.x,axis.y,axis.z = 0.0, 1.0, 0.0
        up_vector = vectors.vector()
        up_vector.x, up_vector.y, up_vector.z = 1.0, 0.0, 0.0
        rmax = window
        region_size = np.array([2.0*rmax,2.0*rmax],order='F',dtype=np.float64)
        if thickness[0] != 0.0:
            thickness_value = group.obj.quantity(thickness[0],thickness[1]).in_units('code_length').d
            distance = 0.5*thickness_value
            far_cut_depth = 0.5*thickness_value
        else:
            distance = rmax
            far_cut_depth = rmax
        enclosing_sphere_p = group.position
        enclosing_sphere_r = rmax
    elif proj.pov == 'z':
        centre = vectors.vector()
        centre.x, centre.y, centre.z = group.position[0].to('code_length')/boxlen, group.position[1].to('code_length')/boxlen, group.position[2].to('code_length')/boxlen
        axis = vectors.vector()
        axis.x,axis.y,axis.z = 0.0, 0.0, 1.0
        up_vector = vectors.vector()
        up_vector.x, up_vector.y, up_vector.z = 1.0, 0.0, 0.0
        rmax = window
        region_size = np.array([2.0*rmax,2.0*rmax],order='F',dtype=np.float64)
        if thickness[0] != 0.0:
            thickness_value = group.obj.quantity(thickness[0],thickness[1]).in_units('code_length').d
            distance = 0.5*thickness_value
            far_cut_depth = 0.5*thickness_value
        else:
            distance = rmax
            far_cut_depth = rmax
        enclosing_sphere_p = group.position
        enclosing_sphere_r = rmax
    elif proj.pov == 'top_midplane':
        axis = vectors.vector()
        norm_L = group.angular_mom['total'].d/np.linalg.norm(group.angular_mom['total'].d)
        los = cartesian_basis['x'] - np.dot(cartesian_basis['x'],norm_L)*norm_L
        los /= np.linalg.norm(los)
        axis.x,axis.y,axis.z = los[0], los[1], los[2]
        up_vector = vectors.vector()
        up_vector.x,up_vector.y,up_vector.z = norm_L[0], norm_L[1], norm_L[2]
        rmax = window
        region_size = np.array([2.0*rmax,2.0*rmax],order='F',dtype=np.float64)
        distance = 0.3*rmax
        far_cut_depth = 0.3*rmax
        centre = vectors.vector()
        im_centre = group.position.d + 0.99*norm_L * rmax
        centre.x, centre.y, centre.z = im_centre[0], im_centre[1], im_centre[2]
        enclosing_sphere_p = im_centre
        enclosing_sphere_r = np.sqrt(max(abs(far_cut_depth),abs(distance))**2 + 2*window**2)
    elif proj.pov == 'bottom_midplane':
        axis = vectors.vector()
        norm_L = -group.angular_mom['total'].d/np.linalg.norm(group.angular_mom['total'].d)
        los = cartesian_basis['x'] - np.dot(cartesian_basis['x'],norm_L)*norm_L
        los /= np.linalg.norm(los)
        axis.x,axis.y,axis.z = los[0], los[1], los[2]
        up_vector = vectors.vector()
        up_vector.x,up_vector.y,up_vector.z = norm_L[0], norm_L[1], norm_L[2]
        rmax = window
        region_size = np.array([2.0*rmax,2.0*rmax],order='F',dtype=np.float64)
        distance = rmax
        far_cut_depth = rmax
        centre = vectors.vector()
        im_centre = group.position.d + norm_L * rmax
        centre.x, centre.y, centre.z = im_centre[0], im_centre[1], im_centre[2]
        enclosing_sphere_p = im_centre
        enclosing_sphere_r = np.sqrt(max(abs(far_cut_depth),abs(distance))**2 + 2*window**2)
    elif proj.pov == 'custom':
        centre = vectors.vector()
        centre.x, centre.y, centre.z = group.position[0], group.position[1], group.position[2]
        axis = vectors.vector()
        norm_L = myaxis/np.linalg.norm(myaxis)
        up = up_axis / np.linalg.norm(up_axis)
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        up_vector = vectors.vector()
        up_vector.x, up_vector.y, up_vector.z = up[0], up[1], up[2]
        rmax = window
        region_size = np.array([2.0*rmax,2.0*rmax],order='F',dtype=np.float64)
        distance = rmax
        far_cut_depth = rmax
        enclosing_sphere_p = group.position
        enclosing_sphere_r = rmax
    elif proj.pov == 'custom_cylinder':
        centre = vectors.vector()
        mycentre = obj.array(mycentre[0],str(mycentre[1])).in_units('code_length')
        centre.x, centre.y, centre.z = mycentre[0].d,mycentre[1].d,mycentre[2].d
        axis = vectors.vector()
        norm_L = myaxis/np.linalg.norm(myaxis)
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        up = cartesian_basis['x'] - np.dot(cartesian_basis['x'],norm_L)*norm_L
        up /= np.linalg.norm(up)
        up_vector = vectors.vector()
        up_vector.x, up_vector.y, up_vector.z = up[0], up[1], up[2]
        region_axis.x,region_axis.y,region_axis.z = norm_L[0], norm_L[1], norm_L[2]
        
        if rmax[1] == 'rvir':
            window = rmax[0]*obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            window = obj.quantity(rmax[0],str(rmax[1])).in_units('code_length')
        rmax = window
        region_size = np.array([2.0*rmax,2.0*rmax],order='F',dtype=np.float64)
        
        if zmin[1] == 'rvir':
            distance = zmin[0]*obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            distance = obj.quantity(zmin[0],str(zmin[1])).in_units('code_length')
        if zmax[1] == 'rvir':
            far_cut_depth = zmax[0]*obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            far_cut_depth = obj.quantity(zmax[0],str(zmax[1])).in_units('code_length')
        enclosing_sphere_p = mycentre + norm_L * max(rmax,0.5*(far_cut_depth-distance))
        enclosing_sphere_r = np.sqrt(max(abs(far_cut_depth),abs(distance))**2 + 2*window**2)
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
        enclosing_sphere_p = group.position
        enclosing_sphere_r = rmax
    # Now create filters if any conditions have been given...
    nfilter = len(filter_conds)
    for i in range(0,nfilter):
        f = init_filter(filter_conds[i],filter_name[i],group)
        proj._get_python_filter(f)
        
    # Do it for hydro first
    filts = []
    for i in range(0,nfilter):
        if isinstance(filter_conds[i],list):
            cond_var = filter_conds[i][0].split('/')[0]
        else:
            cond_var = filter_conds[i].split('/')[0]
        if cond_var in geometrical_variables or cond_var in raw_gas_variables \
            or cond_var in derived_gas_variables or cond_var in gravity_variables:
            f = init_filter(filter_conds[i],filter_name[i],group)
        else:
            # When a filter asks for a variable not existent in the common_variables
            # or the grid_variables dictionaries just ignore it and set it to blank
            f = init_filter('none','none',group)
        filts.append(f)
        
    # Construct substructure regions if we want them out of the projection
    remove_all = False
    if remove_subs == 'all':
        remove_all = True
        remove_subs = True
    if remove_all:
        subs = structure_regions(group, add_substructure=False, add_neighbours=False,
                                    add_intersections=True,position=enclosing_sphere_p,
                                    radius=enclosing_sphere_r)
        nsubs = len(subs)
    elif remove_subs:
        subs = structure_regions(group, add_substructure=True, add_neighbours=False,
                                    tidal_method='BT87_simple')
        nsubs = len(subs)
    else:
        nsubs = 0 
    cam = obs_instruments.init_camera(centre,axis,up_vector,region_size/boxlen,region_axis,bulk,distance/boxlen,
                                      far_cut_depth/boxlen,map_max_size-1,nsubs)
    # Now give filters to camera Fortran type
    if remove_subs and nsubs>0:
        for i in range(0,nsubs):
            cam.subs[i] = subs[i]

    # Update projection details with the camera ones
    proj.width = region_size
    proj.centre = group.position.to('code_length')
    proj.distance = distance/boxlen
    proj.far_cut_depth = far_cut_depth/boxlen
    proj.los_axis = np.array([cam.los_axis.x,cam.los_axis.y,cam.los_axis.z])
    proj.up_vector = np.array([cam.up_vector.x,cam.up_vector.y,cam.up_vector.z])

    # Save the final resolution that will be computed (such as for aspect_ratio>1)
    obs_instruments.get_map_size(cam,proj.resolution)

    # Create projection_handler Fortran derived type for the results of the hydro data projection
    hydro_handler = maps.projection_handler()
    hydro_handler.pov = proj.pov
    hydro_handler.nvars = len(proj.vars['gas'])
    hydro_handler.nwvars = len(proj.weight['gas'])
    hydro_handler.nfilter = nfilter
    maps.allocate_projection_handler(hydro_handler)
    for i in range(0, len(proj.vars['gas'])):
        hydro_handler.varnames.T.view('S128')[i] = proj.vars['gas'][i].ljust(128)
    for i in range(0, len(proj.weight['gas'])):
        hydro_handler.weightvars.T.view('S128')[i] = proj.weight['gas'][i].ljust(128)
    for i in range(0, nfilter):
        hydro_handler.filters[i] = filts[i]
    
    # COMPUTE HYDRO PROJECTION
    if verbose:
        io_ramses.activate_verbose()
    if len(proj.vars['gas']) != 0:
        if verbose: print('Performing hydro projection for '+str(len(proj.vars['gas']))+' variables')
        if obj.use_vardict:
            maps.projection_hydro(group.obj.simulation.fullpath,type_projection,cam,use_neigh,
                                  hydro_handler,int(lmax),int(lmin),nexp_factor,obj.vardict)
        else:
            maps.projection_hydro(group.obj.simulation.fullpath,type_projection,cam,use_neigh,
                                    hydro_handler,int(lmax),int(lmin),nexp_factor)
        # TODO: Weird issue when the direct map array is given.
        data = np.copy(hydro_handler.map)
        proj.data_maps.append(data)
    else:
        del proj.vars['gas']

    # Now settup the filters for particles, making sure
    # we do not include filters inexistent in particle data
    filts = []
    for i in range(0,nfilter):
        if isinstance(filter_conds[i],list):
            cond_var = filter_conds[i][0].split('/')[0]
        else:
            cond_var = filter_conds[i].split('/')[0]
        corrected_var = '_'.join(cond_var.split('_')[1:])
        if corrected_var in geometrical_variables or corrected_var in raw_star_variables \
            or corrected_var in derived_star_variables or corrected_var in raw_dm_variables \
            or corrected_var in derived_dm_variables:
            f = init_filter(filter_conds[i],filter_name[i],group)
        else:
            # When a filter asks for a variable not existent in the common_variables
            # or the particle_variables dictionaries just ignore it and set it to blank
            f = init_filter('none','none',group)
        filts.append(f)
    # Now give filters to camera Fortran type - updating the previous from hydro
    cam = obs_instruments.init_camera(centre,axis,up_vector,region_size/boxlen,region_axis,bulk,distance/boxlen,
                                      far_cut_depth/boxlen,map_max_size,nsubs)
    
    if remove_subs and nsubs>0:
        for i in range(0,nsubs):
            cam.subs[i] = subs[i]
        
    # Create projection_handler Fortran derived type for the results of the hydro data projection
    parts_handler = maps.projection_handler()
    parts_handler.type = proj.pov
    parts_handler.nvars = len(proj.vars['dm'])+len(proj.vars['star'])
    parts_handler.nwvars = len(proj.weight['dm'])+len(proj.weight['star'])
    parts_handler.nfilter = nfilter
    maps.allocate_projection_handler(parts_handler)

    for i in range(0, len(proj.vars['star'])):
        tempstr = 'star/'+proj.vars['star'][i]
        parts_handler.varnames.T.view('S128')[i] = tempstr.ljust(128)
    for i in range(len(proj.vars['star']), len(proj.vars['star'])+len(proj.vars['dm'])):
        tempstr = 'dm/'+proj.vars['dm'][i-len(proj.vars['star'])]
        parts_handler.varnames.T.view('S128')[i] = tempstr.ljust(128)

    for i in range(0, len(proj.weight['star'])):
        tempstr = 'star/'+proj.weight['star'][i]
        parts_handler.weightvars.T.view('S128')[i] = tempstr.ljust(128)
    for i in range(len(proj.weight['star']), len(proj.weight['star'])+len(proj.weight['dm'])):
        tempstr = 'dm/'+proj.weight['dm'][i-len(proj.weight['star'])]
        parts_handler.weightvars.T.view('S128')[i] = tempstr.ljust(128)

    for i in range(0, nfilter):
        parts_handler.filters[i] = filts[i]

    # COMPUTE PARTICLES PROJECTION
    if len(proj.vars['star'])+len(proj.vars['dm']) != 0:
        print('Performing particle projection for '+str(len(proj.vars['star'])+len(proj.vars['dm']))+' variables')
        if tag_file != None and obj.use_part_vardict:
            print('Using the tag file for particles: ',tag_file)
            maps.projection_parts(obj.simulation.fullpath,cam,parts_handler,part_dict=obj.part_vardict,
                                    part_vtypes=part_vartypes,tag_file=tag_file,inverse_tag=inverse_tag)
        elif obj.use_part_vardict:
            maps.projection_parts(obj.simulation.fullpath,cam,parts_handler,part_dict=obj.part_vardict,
                                    part_vtypes=part_vartypes)
        elif tag_file != None:
            print('Using the tag file for particles: ',tag_file)
            maps.projection_parts(obj.simulation.fullpath,cam,parts_handler,tag_file=tag_file,inverse_tag=inverse_tag)
        else:
            maps.projection_parts(obj.simulation.fullpath,cam,parts_handler)
        # TODO: Weird issue when the direct toto array is given.
        data = np.copy(parts_handler.toto)
        if len(proj.vars['star']) == 0:
            data_dm = data.reshape(nfilter,len(proj.vars['dm']),data.shape[2],data.shape[3])
            proj.data_maps.append(data_dm)
        elif len(proj.vars['dm']) == 0:
            data_star = data[:,:len(proj.vars['star'])].reshape(nfilter,len(proj.vars['star']),data.shape[2],data.shape[3])
            proj.data_maps.append(data_star)
        else:
            data_star = data[:,:len(proj.vars['star'])].reshape(nfilter,len(proj.vars['star']),data.shape[2],data.shape[3])
            proj.data_maps.append(data_star)
            data_dm = data[:,len(proj.vars['star']):].reshape(nfilter,len(proj.vars['dm']),data.shape[2],data.shape[3])
            proj.data_maps.append(data_dm)
    else:
        if 'star' in proj.vars.keys():
            del proj.vars['star']
        if 'dm' in proj.vars.keys():
            del proj.vars['dm']
    if 'star' in proj.vars.keys():
        if len(proj.vars['star']) == 0:
            del proj.vars['star']
    if 'dm' in proj.vars.keys():
        if len(proj.vars['dm']) == 0:
            del proj.vars['dm']

    return proj


def plot_single_galaxy_projection(proj_FITS,fields,logscale=True,scalebar=(3,'kpc'),redshift=True,returnfig=False,
                                    pov='x',centers=[],radii=[],circle_keys=[],filter_name='none',smooth=False,
                                    type_scale='galaxy'):
    """This function uses the projection information in a FITS file following the 
        OZY format and plots it following the OZY standards."""
    
    # Make required imports
    from ozy.utils import tidal_radius,invert_tick_colours
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from matplotlib.colors import LogNorm,SymLogNorm
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import seaborn as sns
    from unyt import unyt_quantity,unyt_array
    import scipy as sp                                                                                                                                                                                     
    import scipy.ndimage                                                                                                                                                                                   
    sns.set(style="dark")
    plt.rcParams["axes.axisbelow"] = False
    plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "Computer Modern Roman",
    })

    if type_scale != '':
        type_scale = '_'+type_scale

    # First,check that FITS file actually exists
    if not os.path.exists(proj_FITS):
        raise ImportError('File not found. Please check!: '+str(proj_FITS))


    # Load FITS file
    hdul = fits.open(proj_FITS)
    hdul_fields = [h.header['btype'] for h in hdul]

    # Add filter identification to the field
    hdul_filtername = ['COND_0' in h.header for h in hdul]
    if any(hdul_filtername):
        for f in range(0,len(fields)):
            fields[f] = fields[f] + '/' + filter_name
    else:
        print('WARNING: This FITS file does not have filter information, so anything you added will be ignored!')

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

    # Construct camera
    length_unit = unyt_quantity(1.0,hdul[0].header['CUNIT1'])
    width_x =  hdul[0].header['CDELT1']*hdul[0].header['NAXIS1'] * length_unit
    width_y =  hdul[0].header['CDELT2']*hdul[0].header['NAXIS2'] * length_unit
    ex = [-0.5*width_x,0.5*width_x,-0.5*width_y,0.5*width_y]
    los_axis,up_axis,centre,velocity = vectors.vector(),vectors.vector(),vectors.vector(),vectors.vector()
    los_axis.x,los_axis.y,los_axis.z = hdul[0].header['LOS_X'],hdul[0].header['LOS_Y'],hdul[0].header['LOS_Z']
    up_axis.x,up_axis.y,up_axis.z = hdul[0].header['UP_X'],hdul[0].header['UP_Y'],hdul[0].header['UP_Z']
    centre.x,centre.y,centre.z = hdul[0].header['CENTRE_X'],hdul[0].header['CENTRE_Y'],hdul[0].header['CENTRE_Z']
    cam = obs_instruments.init_camera(centre,los_axis,up_axis,width_x.d,los_axis,velocity,
                                      width_x.d,width_x.d,hdul[0].header['NAXIS1'],
                                      1,len(centers))
    stellar = False
    for i in range(0, ax.shape[0]):
        for j in range(0, ax.shape[1]):
            ivar = i*ax.shape[1] + j
            if ivar >= len(fields):
                # Clear that extra panel
                ax[i,j].get_xaxis().set_visible(False)
                ax[i,j].get_yaxis().set_visible(False)
                break
            ax[i,j].set_xlim([-0.5*width_x,0.5*width_x])
            ax[i,j].set_ylim([-0.5*width_y,0.5*width_y])
            ax[i,j].axes.xaxis.set_visible(False)
            ax[i,j].axes.yaxis.set_visible(False)
            ax[i,j].axis('off')
            h = [k for k in range(0,len(hdul)) if hdul[k].header['btype']==fields[ivar]][0]

            if fields[ivar].split('/')[0] == 'star' or fields[ivar].split('/')[0] == 'dm':
                plotting_def = get_plotting_def(fields[ivar].split('/')[0]+'_'+fields[ivar].split('/')[1])
                stellar = True
                full_varname = fields[ivar].split('/')[0]+'_'+fields[ivar].split('/')[1]
            else:
                plotting_def = get_plotting_def(fields[ivar].split('/')[1])
                full_varname = fields[ivar].split('/')[1]
            print(fields[ivar],np.nanmin(hdul[h].data.T),np.nanmax(hdul[h].data.T))
            if smooth:
                # We smooth with a Gaussian kernel with x_stddev=1 (and y_stddev=1)
                # It is a 9x9 array
                from astropy.convolution import Gaussian2DKernel, convolve
                kernel = Gaussian2DKernel(x_stddev=2.5)
                # astropy's convolution replaces the NaN pixels with a kernel-weighted
                # interpolation from their neighbors
                cImage = convolve(hdul[h].data.T, kernel)
                # sigma=3
                # cImage = sp.ndimage.filters.gaussian_filter(hdul[h].data.T, sigma, mode='constant')
            else:
                cImage = hdul[h].data.T
            if logscale and not plotting_def['symlog']:
                
                plot = ax[i,j].imshow(np.log10(cImage), cmap=plotting_def['cmap'],
                                origin='upper',vmin=np.log10(plotting_def['vmin'+type_scale]),
                                vmax=np.log10(plotting_def['vmax'+type_scale]),extent=ex,
                                interpolation='nearest')
            elif logscale and plotting_def['symlog']:
                if (filter_name != 'outflow' and filter_name != 'inflow'):
                    plot = ax[i,j].imshow(cImage, cmap=plotting_def['cmap'],
                                    origin='upper',norm=SymLogNorm(linthresh=plotting_def['linthresh'], linscale=1,
                                    vmin=plotting_def['vmin'+type_scale], vmax=plotting_def['vmax'+type_scale]),
                                    extent=ex,
                                    interpolation='none')
                elif fields[ivar].split('/')[1] == 'v_sphere_r':
                    if filter_name == 'inflow' and smooth:                                                                                                                                                                                          
                        cImage = sp.ndimage.filters.gaussian_filter(-hdul[h].data.T, sigma, mode='constant')
                    plot = ax[i,j].imshow(np.log10(cImage), cmap=plotting_def['cmap_'+filter_name],
                                    origin='upper',vmin=np.log10(plotting_def['vmin_'+filter_name]),
                                    vmax=np.log10(plotting_def['vmax_'+filter_name]),extent=ex,
                                    interpolation='none')
                else:
                    plot = ax[i,j].imshow(cImage, cmap=plotting_def['cmap'],
                                    origin='upper',norm=SymLogNorm(linthresh=plotting_def['linthresh'], linscale=1,
                                    vmin=plotting_def['vmin'+type_scale], vmax=plotting_def['vmax'+type_scale]),
                                    extent=ex,
                                    interpolation='none')
            else:
                plot = ax[i,j].imshow(hdul[h].data.T, cmap=plotting_def['cmap'],
                                origin='upper',extent=ex,interpolation='none',
                                vmin=plotting_def['vmin'+type_scale],vmax=plotting_def['vmax'+type_scale])

            cbaxes = inset_axes(ax[i,j], width="80%", height="6%", loc='lower center')
            cbar = fig.colorbar(plot, cax=cbaxes, orientation='horizontal')
            if logscale and not plotting_def['symlog']:
                cbar.set_label(plotting_def['label_log'],color=plotting_def['text_over'],fontsize=12,labelpad=-27, y=0.85,weight='bold')
            else:
                cbar.set_label(plotting_def['label'],color=plotting_def['text_over'],fontsize=12,labelpad=-27, y=0.85)
            cbar.ax.tick_params(axis='x', pad=-8.05, labelsize=11,labelcolor=plotting_def['text_over'])
            cbar.ax.tick_params(length=0,width=0)
            cbar.outline.set_linewidth(1)
            cbar.ax.set_axisbelow(False)
            invert_tick_colours(cbar.ax,full_varname,type_scale)

            if redshift and i==0 and j==0:
                ax[i,j].text(0.05, 0.90, r'$z = {z:.2f}$'.format(z=hdul[h].header['redshift']), # Redshift
                                    verticalalignment='bottom', horizontalalignment='left',
                                    transform=ax[i,j].transAxes,
                                    color=plotting_def['text_over'], fontsize=10,fontweight='bold')

            fontprops = fm.FontProperties(size=10,weight='bold')
            if scalebar != None:
                if scalebar and i==0 and j==0:
                    sl = int(unyt_quantity(scalebar[0],scalebar[1]).in_units(hdul[0].header['CUNIT1']).d)
                    slb = AnchoredSizeBar(ax[i,j].transData,
                                            sl, '%s %s'%(str(int(sl)),hdul[0].header['CUNIT1']), 
                                            'upper right', 
                                            pad=0.1,
                                            color=plotting_def['text_over'],
                                            frameon=False,
                                            size_vertical=0.1,
                                            fontproperties=fontprops)
                    ax[i,j].add_artist(slb)

            if len(centers) != 0 and len(radii) != 0 and len(circle_keys) != 0:
                # This expects for each point:
                # 1. a center given as a OZY.array
                # 2. a radius given as an OZY.quantity
                # 3. a dictionary of the plotting settings of the circle,
                #    including the "edgecolor" and the "linestyle"
                for c in range(0, len(centers)):
                    centrecircle = centers[c].in_units(hdul[0].header['CUNIT1']).d
                    dist = centrecircle - np.array([centre.x,centre.y,centre.z])
                    dist = np.linalg.norm(dist) - radii[c].in_units(hdul[0].header['CUNIT1']).d
                    if dist <= (width_x/2):
                        obs_instruments.project_points(cam,1,centrecircle)
                        r = radii[c].in_units(hdul[0].header['CUNIT1']).d
                        circle_settings = circle_keys[c]
                        circle = plt.Circle(centrecircle,r,fill=False,
                                            edgecolor=circle_settings['edgecolor'],
                                            linestyle=circle_settings['linestyle'],
                                            linewidth=circle_settings['linewidth'])
                        ax[i,j].add_patch(circle)
    hdul.close()
    fig.subplots_adjust(hspace=0,wspace=0,left=0,right=1, bottom=0, top=1)
    if returnfig:
        return fig
    else:
        if filter_name != 'none':
            fig.savefig(proj_FITS.split('.')[0]+'_'+filter_name+'.png',format='png',dpi=300)
        else:
            fig.savefig(proj_FITS.split('.')[0]+'.png',format='png',dpi=300)

def plot_fe(faceon_fits,edgeon_fits,fields,logscale=True,scalebar=(3,'kpc'),
                 redshift=True,returnfig=False,pov='x',type_scale='galaxy',
                 smooth=False,
                 labels=[],centers=[],radii=[],circle_keys=[]):
    """This function uses the projection information in a set of FITS files and combines in a single
        figure following the OZY standards."""
    
    # Make required imports
    from ozy.utils import tidal_radius,invert_tick_colours
    from ozy.plot_settings import circle_dictionary
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from matplotlib.colors import LogNorm,SymLogNorm
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import seaborn as sns
    from unyt import unyt_quantity,unyt_array
    import scipy as sp                                                                                                                                                                                     
    import scipy.ndimage
    sns.set(style="dark")
    plt.rcParams["axes.axisbelow"] = False
    plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "Computer Modern Roman",
    })
    if type_scale != '':
        type_scale = '_'+type_scale

    # First,check that the FITS files actually exist
    if not os.path.exists(faceon_fits):
        raise ImportError('File not found. Please check!: '+str(faceon_fits))
    if not os.path.exists(edgeon_fits):
        raise ImportError('File not found. Please check!: '+str(edgeon_fits))

    # Since everything is fine, we begin plotting…
    ncol = int(len(fields))
    nrow = 2
    height_ratios = [2,1]
    figsize = plt.figaspect(3 / (2 * float(ncol)))
    fig = plt.figure(figsize=figsize, facecolor='k', edgecolor='k')
    plot_grid = fig.add_gridspec(nrow, ncol, wspace=0, hspace=0,left=0,right=1,
                                    bottom=0, top=1, height_ratios=height_ratios)
    ax = []
    for i in range(0,nrow):
        ax.append([])
        for j in range(0,ncol):
            ax[i].append(fig.add_subplot(plot_grid[i,j]))
    ax = np.asarray(ax)

    hdul = [fits.open(faceon_fits), fits.open(edgeon_fits)]
    
    for j in range(0, 2):
        for i in range(0, len(fields)):
            # Construct camera
            length_unit = unyt_quantity(1.0,hdul[j][0].header['CUNIT1'])
            width_x =  hdul[j][0].header['CDELT1']*hdul[j][0].header['NAXIS1'] * length_unit
            width_y =  hdul[j][0].header['CDELT2']*hdul[j][0].header['NAXIS2'] * length_unit
            ex = [-0.5*width_x,0.5*width_x,-0.5*width_y,0.5*width_y]
            los_axis,up_axis,centre,velocity = vectors.vector(),vectors.vector(),vectors.vector(),vectors.vector()
            los_axis.x,los_axis.y,los_axis.z = hdul[j][0].header['LOS_X'],hdul[j][0].header['LOS_Y'],hdul[j][0].header['LOS_Z']
            up_axis.x,up_axis.y,up_axis.z = hdul[j][0].header['UP_X'],hdul[j][0].header['UP_Y'],hdul[j][0].header['UP_Z']
            centre.x,centre.y,centre.z = hdul[j][0].header['CENTRE_X'],hdul[j][0].header['CENTRE_Y'],hdul[j][0].header['CENTRE_Z']
            cam = obs_instruments.init_camera(centre,los_axis,up_axis,width_x.d,los_axis,velocity,
                                            width_x.d,width_x.d,hdul[j][0].header['NAXIS1'],
                                            1,len(centers))
            if j%2 == 0:
                ax[j,i].set_xlim([-0.5*width_x,0.5*width_x])
                ax[j,i].set_ylim([-0.5*width_y,0.5*width_y])
                ex = np.asarray([-0.5*width_x,0.5*width_x,-0.5*width_y,0.5*width_y])
            else:
                ax[j,i].set_xlim([-0.5*width_x,0.5*width_x])
                ax[j,i].set_ylim([-0.25*width_y,0.25*width_y])
                ex = np.asarray([-0.5*width_x,0.5*width_x,-0.5*width_y,0.5*width_y])
            ax[j,i].axes.xaxis.set_visible(False)
            ax[j,i].axes.yaxis.set_visible(False)
            ax[j,i].axis('off')
            if fields[i].split('/')[1] != 'densityandv_sphere_r':
                h = [k for k in range(0,len(hdul[j])) if hdul[j][k].header['btype']==fields[i]][0]

            if fields[i].split('/')[0] == 'star' or fields[i].split('/')[0] == 'dm':
                plotting_def = get_plotting_def(fields[i].split('/')[0]+'_'+fields[i].split('/')[1])
                full_varname = fields[i].split('/')[0]+'_'+fields[i].split('/')[1]
            elif fields[i].split('/')[1] != 'densityandv_sphere_r':
                plotting_def = get_plotting_def(fields[i].split('/')[1])
                full_varname = fields[i].split('/')[1]
            if smooth:
                # We smooth with a Gaussian kernel with x_stddev=1 (and y_stddev=1)
                # It is a 9x9 array
                from astropy.convolution import Gaussian2DKernel, convolve
                kernel = Gaussian2DKernel(x_stddev=1.5)
                # astropy's convolution replaces the NaN pixels with a kernel-weighted
                # interpolation from their neighbors
                cImage = convolve(hdul[j][h].data.T, kernel)
                # sigma=3
                # cImage = sp.ndimage.filters.gaussian_filter(hdul[j][h].data.T, sigma, mode='constant')
            else:
                cImage = hdul[j][h].data.T
                
                        
            if logscale and fields[i].split('/')[1] != 'v_sphere_r':
                plot = ax[j,i].imshow(np.log10(cImage), cmap=plotting_def['cmap'],
                                origin='upper',vmin=np.log10(plotting_def['vmin'+type_scale]),
                                vmax=np.log10(plotting_def['vmax'+type_scale]),extent=ex,
                                interpolation='nearest')
            elif logscale and fields[i].split('/')[1] == 'v_sphere_r':
                plot = ax[j,i].imshow(cImage, cmap=plotting_def['cmap'],
                                origin='upper',norm=SymLogNorm(linthresh=10, linscale=1,vmin=plotting_def['vmin'], vmax=plotting_def['vmax']),
                                extent=ex,
                                interpolation='nearest')
            
            else:
                plot = ax[j,i].imshow(cImage, cmap=plotting_def['cmap'],
                                origin='upper',extent=ex,interpolation='nearest',
                                vmin=plotting_def['vmin'],vmax=plotting_def['vmax'])
            if j%2 == 0:
                cbaxes = inset_axes(ax[j,i], width="80%", height="6%", loc='lower center')
                cbar = fig.colorbar(plot, cax=cbaxes, orientation='horizontal')
                if logscale and not plotting_def['symlog']:
                    cbar.set_label(plotting_def['label_log'],color=plotting_def['text_over'],fontsize=14,labelpad=-30, y=0.85,weight='bold')
                else:
                    cbar.set_label(plotting_def['label'],color=plotting_def['text_over'],fontsize=14,labelpad=-30, y=0.85)
                cbar.ax.tick_params(axis='x', pad=-8.05, labelsize=11,labelcolor=plotting_def['text_over'])
                cbar.ax.tick_params(length=0,width=0)
                cbar.outline.set_linewidth(1)
                cbar.ax.set_axisbelow(False)
                invert_tick_colours(cbar.ax,full_varname,type_scale)
            if redshift and j%2==0 and i==0:
                ax[j,i].text(0.05, 0.88, r'$z = {z:.2f}$'.format(z=hdul[j][0].header['redshift']), # Redshift
                                    verticalalignment='bottom', horizontalalignment='left',
                                    transform=ax[j,i].transAxes,
                                    color=plotting_def['text_over'], fontsize=15,fontweight='bold')

            fontprops = fm.FontProperties(size=15,weight='bold')
            if scalebar != None:
                if scalebar and j%2==0 and i==0:
                    sl = int(unyt_quantity(scalebar[0],scalebar[1]).in_units(hdul[j][0].header['CUNIT1']).d)
                    slb = AnchoredSizeBar(ax[j,i].transData,
                                                sl, '%s %s'%(str(int(sl)),hdul[j][0].header['CUNIT1']),
                                                'upper right',
                                                pad=0.1,
                                                color=plotting_def['text_over'],
                                                frameon=False,
                                                size_vertical=0.1,
                                                fontproperties=fontprops)
                    ax[j,i].add_artist(slb)
            if len(labels) !=0 and j%2==0 and i==ax.shape[1]-1:
                l = labels[int(i/2)][0] + '\n' + labels[int(i/2)][1] + '\n' + labels[int(i/2)][2]
                ax[j,i].text(0.15, 0.7, l,
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax[j,i].transAxes,
                                color=plotting_def['text_over'], fontsize=10,fontweight='bold',
                                bbox=dict(facecolor='none', edgecolor='white'))
            if len(centers) != 0 and len(radii) != 0 and len(circle_keys) != 0 and i%2 == 0:
                # This expects for each point:
                # 1. a center given as a OZY.array
                # 2. a radius given as an OZY.quantity
                # 3. a dictionary of the plotting settings of the circle,
                #    including the "edgecolor" and the "linestyle"
                centrecircle = centers[int(i/2)].in_units(hdul[j][0].header['CUNIT1']).d
                dist = centrecircle - np.array([centre.x,centre.y,centre.z])
                dist = np.linalg.norm(dist) - radii[int(i/2)].in_units(hdul[j][0].header['CUNIT1']).d
                if dist <= (width_x/2):
                    obs_instruments.project_points(cam,1,centrecircle)
                    r = radii[int(i/2)].in_units(hdul[j][0].header['CUNIT1']).d
                    circle_settings = circle_keys[int(i/2)]
                    circle = plt.Circle(centrecircle,r,fill=False,
                                        edgecolor=circle_settings['edgecolor'],
                                        linestyle=circle_settings['linestyle'],
                                        linewidth=circle_settings['linewidth'],
                                        alpha=1.0)
                    ax[j,i].add_patch(circle)

    fig.subplots_adjust(hspace=0,wspace=0,left=0,right=1, bottom=0, top=1)
    if returnfig:
        return fig
    else:
        fig.savefig('plot_facevsedge.png',format='png',dpi=300)

def plot_comp_fe(faceon_fits,edgeon_fits,fields,logscale=True,scalebar=(3,'kpc'),
                 redshift=True,returnfig=False,pov='x',type_scale='galaxy',
                 skirt_images=None,smooth=False,
                 labels=[],centers=[],radii=[],circle_keys=[]):
    """This function uses the projection information in a set of FITS files and combines in a single
        figure following the OZY standards."""
    
    # Make required imports
    from ozy.utils import tidal_radius,invert_tick_colours
    from ozy.variables_settings import circle_dictionary
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from matplotlib.colors import LogNorm,SymLogNorm
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import seaborn as sns
    from unyt import unyt_quantity,unyt_array
    import scipy as sp                                                                                                                                                                                     
    import scipy.ndimage
    sns.set(style="dark")
    plt.rcParams["axes.axisbelow"] = False
    plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "Computer Modern Roman",
    })
    if type_scale != '':
        type_scale = '_'+type_scale

    # Check if we have the same number of faceon and edgeon files
    if len(faceon_fits) != len(edgeon_fits):
        raise ImportError('The number of faceon to edgeon FITS is not the same: %i to %i'%(len(faceon_fits),len(edgeon_fits)))

    # First,check that the FITS files actually exist
    for f in faceon_fits:
        if not os.path.exists(f):
            raise ImportError('File not found. Please check!: '+str(f))
    for f in edgeon_fits:
        if not os.path.exists(f):
            raise ImportError('File not found. Please check!: '+str(f))

    # Since everything is fine, we begin plotting…
    ncol = int(len(fields))
    nrow = int(len(faceon_fits)*2)
    height_ratios = []
    for i in range(0, len(faceon_fits)):
        height_ratios.append(2)
        height_ratios.append(1)
    figsize = plt.figaspect((3 * float(len(faceon_fits))) / (2 * float(ncol)))
    fig = plt.figure(figsize=2*figsize, facecolor='k', edgecolor='k')
    plot_grid = fig.add_gridspec(nrow, ncol, wspace=0, hspace=0,left=0,right=1,
                                    bottom=0, top=1, height_ratios=height_ratios)
    ax = []
    for i in range(0,nrow):
        ax.append([])
        for j in range(0,ncol):
            ax[i].append(fig.add_subplot(plot_grid[i,j]))
            # if i==0 and j==0:
            #     axmain = fig.add_subplot(plot_grid[i,j])
            #     ax[i].append(axmain)
            # else:
            #     ax[i].append(fig.add_subplot(plot_grid[i,j], sharex=axmain, sharey=axmain))
    ax = np.asarray(ax)
    
    for i in range(0, ax.shape[0]):
        if i%2 == 0:
            myfits = faceon_fits[int(i/2)]
        else:
            myfits = edgeon_fits[int(i/2)]
        hdul = fits.open(myfits)
        hdul_fields = [h.header['btype'] for h in hdul]
        # Construct camera
        length_unit = unyt_quantity(1.0,hdul[0].header['CUNIT1'])
        width_x =  hdul[0].header['CDELT1']*hdul[0].header['NAXIS1'] * length_unit
        width_y =  hdul[0].header['CDELT2']*hdul[0].header['NAXIS2'] * length_unit
        ex = [-0.5*width_x,0.5*width_x,-0.5*width_y,0.5*width_y]
        los_axis,up_axis,centre,velocity = vectors.vector(),vectors.vector(),vectors.vector(),vectors.vector()
        los_axis.x,los_axis.y,los_axis.z = hdul[0].header['LOS_X'],hdul[0].header['LOS_Y'],hdul[0].header['LOS_Z']
        up_axis.x,up_axis.y,up_axis.z = hdul[0].header['UP_X'],hdul[0].header['UP_Y'],hdul[0].header['UP_Z']
        centre.x,centre.y,centre.z = hdul[0].header['CENTRE_X'],hdul[0].header['CENTRE_Y'],hdul[0].header['CENTRE_Z']
        cam = obs_instruments.init_camera(centre,los_axis,up_axis,width_x.d,los_axis,velocity,
                                        width_x.d,width_x.d,hdul[0].header['NAXIS1'],
                                        1,len(centers))
        
        for j in range(0, ax.shape[1]):
            ivar = j
            have_skirt = False
            if fields[ivar] not in hdul_fields and \
                fields[ivar].split('/')[1] != 'densityandv_sphere_r'and \
                fields[ivar].split('/')[1] != 'skirt':
                # Clear that extra panel
                ax[i,j].get_xaxis().set_visible(False)
                ax[i,j].get_yaxis().set_visible(False)
                continue 
            if fields[ivar].split('/')[1] == 'skirt': have_skirt = True
            if len(skirt_images) != len(faceon_fits):
                print('You need to provide the same number of SKIRT projections than simulations!')
                exit
            if i%2 == 0:
                ax[i,j].set_xlim([-0.5*width_x,0.5*width_x])
                ax[i,j].set_ylim([-0.5*width_y,0.5*width_y])
                ex = np.asarray([-0.5*width_x,0.5*width_x,-0.5*width_y,0.5*width_y])
            else:
                ax[i,j].set_xlim([-0.5*width_x,0.5*width_x])
                ax[i,j].set_ylim([-0.25*width_y,0.25*width_y])
                ex = np.asarray([-0.5*width_x,0.5*width_x,-0.5*width_y,0.5*width_y])
            ax[i,j].axes.xaxis.set_visible(False)
            ax[i,j].axes.yaxis.set_visible(False)
            ax[i,j].axis('off')
            if not have_skirt:
                if fields[ivar].split('/')[1] != 'densityandv_sphere_r':
                    h = [k for k in range(0,len(hdul)) if hdul[k].header['btype']==fields[ivar]][0]

                if fields[ivar].split('/')[0] == 'star' or fields[ivar].split('/')[0] == 'dm':
                    plotting_def = get_plotting_def(fields[ivar].split('/')[0]+'_'+fields[ivar].split('/')[1])
                    full_varname = fields[ivar].split('/')[0]+'_'+fields[ivar].split('/')[1]
                elif fields[ivar].split('/')[1] != 'densityandv_sphere_r':
                    plotting_def = get_plotting_def(fields[ivar].split('/')[1])
                    full_varname = fields[ivar].split('/')[1]
                if smooth:
                    # We smooth with a Gaussian kernel with x_stddev=1 (and y_stddev=1)
                    # It is a 9x9 array
                    from astropy.convolution import Gaussian2DKernel, convolve
                    kernel = Gaussian2DKernel(x_stddev=1.5)
                    # astropy's convolution replaces the NaN pixels with a kernel-weighted
                    # interpolation from their neighbors
                    cImage = convolve(hdul[h].data.T, kernel)
                    # sigma=3
                    # cImage = sp.ndimage.filters.gaussian_filter(hdul[h].data.T, sigma, mode='constant')
                else:
                    cImage = hdul[h].data.T
                if fields[ivar].split('/')[1] == 'densityandv_sphere_r':
                    plotting_def = get_plotting_def('density')
                    full_varname = 'density'
                    if i%2 == 0:
                        plot = ax[i,j].imshow(np.log10(hdul['gas/density'].data.T), cmap=plotting_def['cmap'],
                                    origin='upper',vmin=np.log10(plotting_def['vmin'+type_scale]),
                                    vmax=np.log10(plotting_def['vmax'+type_scale]),extent=ex,
                                    interpolation='nearest')
                    else:
                        kernel = np.outer(signal.gaussian(100,1),signal.gaussian(100,1))
                        density = np.log10(hdul['gas/density'].data.T)
                        v_sphere = hdul['gas/v_sphere_r'].data.T
                        outflow = np.copy(density)
                        inflow = np.copy(density)
                        outflow[v_sphere <= -1.0] = np.log10(plotting_def['vmin'+type_scale])
                        inflow[v_sphere > 1.0] = np.log10(plotting_def['vmin'+type_scale])
                        outflow = (outflow-np.log10(plotting_def['vmin'+type_scale]))/(np.log10(plotting_def['vmax'+type_scale]) - np.log10(plotting_def['vmin'+type_scale]))
                        inflow = (inflow-np.log10(plotting_def['vmin'+type_scale]))/(np.log10(plotting_def['vmax'+type_scale]) - np.log10(plotting_def['vmin'+type_scale]))
                        print('Density maps coloured with radial velocity')
                        nx,ny = outflow.shape
                        rgbArray = np.zeros((nx,ny,3), 'uint8')
                        rgbArray[:,:,0] = outflow*255*0.8
                        rgbArray[:,:,1] = 0.1*255*0.8
                        rgbArray[:,:,2] = inflow*255*0.8
                        plot = ax[i,j].imshow(np.zeros((nx,ny,3), 'uint8'), cmap='gray',
                                        vmin=plotting_def['vmin'+type_scale],vmax=plotting_def['vmax'+type_scale])
                        ax[i,j].imshow(v_sphere,origin='upper',extent=ex,interpolation='lanczos')
                        
                elif logscale and fields[ivar].split('/')[1] != 'v_sphere_r':
                    plot = ax[i,j].imshow(np.log10(cImage), cmap=plotting_def['cmap'],
                                    origin='upper',vmin=np.log10(plotting_def['vmin'+type_scale]),
                                    vmax=np.log10(plotting_def['vmax'+type_scale]),extent=ex,
                                    interpolation='nearest')
                elif logscale and fields[ivar].split('/')[1] == 'v_sphere_r':
                    plot = ax[i,j].imshow(cImage, cmap=plotting_def['cmap'],
                                    origin='upper',norm=SymLogNorm(linthresh=10, linscale=1,vmin=plotting_def['vmin'], vmax=plotting_def['vmax']),
                                    extent=ex,
                                    interpolation='nearest')
                
                else:
                    plot = ax[i,j].imshow(cImage, cmap=plotting_def['cmap'],
                                    origin='upper',extent=ex,interpolation='nearest',
                                    vmin=plotting_def['vmin'],vmax=plotting_def['vmax'])
                if i%2 == 0:
                    cbaxes = inset_axes(ax[i,j], width="80%", height="6%", loc='lower center')
                    cbar = fig.colorbar(plot, cax=cbaxes, orientation='horizontal')
                    if logscale and not plotting_def['symlog']:
                        cbar.set_label(plotting_def['label_log'],color=plotting_def['text_over'],fontsize=14,labelpad=-30, y=0.85,weight='bold')
                    else:
                        cbar.set_label(plotting_def['label'],color=plotting_def['text_over'],fontsize=14,labelpad=-30, y=0.85)
                    cbar.ax.tick_params(axis='x', pad=-8.05, labelsize=11,labelcolor=plotting_def['text_over'])
                    cbar.ax.tick_params(length=0,width=0)
                    cbar.outline.set_linewidth(1)
                    cbar.ax.set_axisbelow(False)
                    invert_tick_colours(cbar.ax,full_varname,type_scale)
                if redshift and i%2==0 and j==0:
                    ax[i,j].text(0.05, 0.88, r'$z = {z:.2f}$'.format(z=hdul[0].header['redshift']), # Redshift
                                        verticalalignment='bottom', horizontalalignment='left',
                                        transform=ax[i,j].transAxes,
                                        color=plotting_def['text_over'], fontsize=15,fontweight='bold')

                fontprops = fm.FontProperties(size=15,weight='bold')
                if scalebar != None:
                    if scalebar and i%2==0 and j==0:
                        sl = int(unyt_quantity(scalebar[0],scalebar[1]).in_units(hdul[0].header['CUNIT1']).d)
                        slb = AnchoredSizeBar(ax[i,j].transData,
                                                    sl, '%s %s'%(str(int(sl)),hdul[0].header['CUNIT1']),
                                                    'upper right',
                                                    pad=0.1,
                                                    color=plotting_def['text_over'],
                                                    frameon=False,
                                                    size_vertical=0.1,
                                                    fontproperties=fontprops)
                        ax[i,j].add_artist(slb)
                if len(labels) !=0 and i%2==0 and j==ax.shape[1]-1:
                    l = labels[int(i/2)][0] + '\n' + labels[int(i/2)][1] + '\n' + labels[int(i/2)][2]
                    ax[i,j].text(0.15, 0.7, l,
                                    verticalalignment='bottom', horizontalalignment='left',
                                    transform=ax[i,j].transAxes,
                                    color=plotting_def['text_over'], fontsize=10,fontweight='bold',
                                    bbox=dict(facecolor='none', edgecolor='white'))
                if len(centers) != 0 and len(radii) != 0 and len(circle_keys) != 0 and i%2 == 0:
                    # This expects for each point:
                    # 1. a center given as a OZY.array
                    # 2. a radius given as an OZY.quantity
                    # 3. a dictionary of the plotting settings of the circle,
                    #    including the "edgecolor" and the "linestyle"
                    centrecircle = centers[int(i/2)].in_units(hdul[0].header['CUNIT1']).d
                    dist = centrecircle - np.array([centre.x,centre.y,centre.z])
                    dist = np.linalg.norm(dist) - radii[int(i/2)].in_units(hdul[0].header['CUNIT1']).d
                    if dist <= (width_x/2):
                        obs_instruments.project_points(cam,1,centrecircle)
                        r = radii[int(i/2)].in_units(hdul[0].header['CUNIT1']).d
                        circle_settings = circle_keys[int(i/2)]
                        circle = plt.Circle(centrecircle,r,fill=False,
                                            edgecolor=circle_settings['edgecolor'],
                                            linestyle=circle_settings['linestyle'],
                                            linewidth=circle_settings['linewidth'],
                                            alpha=1.0)
                        ax[i,j].add_patch(circle)
                    
            else:
                if i%2 == 0:
                    ski_im = plt.imread(skirt_images[int(i/2)][0])
                    ex = np.asarray([-0.5*width_x,0.5*width_x,-0.5*width_y,0.5*width_y])
                    ax[i,j].imshow(ski_im,origin='upper',extent=-ex)
                else:
                    ski_im = plt.imread(skirt_images[int(i/2)][1])
                    ex = np.asarray([-0.5*width_x,0.5*width_x,-0.5*width_y,0.5*width_y])
                    ax[i,j].imshow(ski_im,origin='upper',extent=ex)

                
                if redshift and i%2==0 and j==0:
                    ax[i,j].text(0.05, 0.90, r'$z = {z:.2f}$'.format(z=hdul[0].header['redshift']), # Redshift
                                        verticalalignment='bottom', horizontalalignment='left',
                                        transform=ax[i,j].transAxes,
                                        color='w', fontsize=15,fontweight='bold')

                fontprops = fm.FontProperties(size=15,weight='bold')
                if scalebar != None:
                    if scalebar and i%2==0:
                        sl = int(unyt_quantity(scalebar[0],scalebar[1]).in_units(hdul[0].header['CUNIT1']).d)
                        slb = AnchoredSizeBar(ax[i,j].transData,
                                                    2*sl, '%s %s'%(str(int(sl)),hdul[0].header['CUNIT1']),
                                                    'lower right',
                                                    pad=0.1,
                                                    color='white',
                                                    frameon=False,
                                                    size_vertical=0.1,
                                                    fontproperties=fontprops)
                        ax[i,j].add_artist(slb)
                if len(labels) !=0 and i%2==0 and j==ax.shape[1]-1:
                    l = labels[int(i/2)][0] + '\n' + labels[int(i/2)][1] + '\n' + labels[int(i/2)][2]
                    ax[i,j].text(0.15, 0.7, l,
                                    verticalalignment='bottom', horizontalalignment='left',
                                    transform=ax[i,j].transAxes,
                                    color='w', fontsize=10,fontweight='bold')
                if len(centers) != 0 and len(radii) != 0 and len(circle_keys) != 0 and i%2 == 0:
                    # This expects for each point:
                    # 1. a center given as a OZY.array
                    # 2. a radius given as an OZY.quantity
                    # 3. a dictionary of the plotting settings of the circle,
                    #    including the "edgecolor" and the "linestyle"
                    centrecircle = centers[int(i/2)].in_units(hdul[0].header['CUNIT1']).d
                    dist = centrecircle - np.array([centre.x,centre.y,centre.z])
                    dist = np.linalg.norm(dist) - 2*radii[int(i/2)].in_units(hdul[0].header['CUNIT1']).d
                    if dist <= (width_x/2):
                        obs_instruments.project_points(cam,1,centrecircle)
                        r = 2*radii[int(i/2)].in_units(hdul[0].header['CUNIT1']).d
                        circle_settings = circle_keys[int(i/2)]
                        circle = plt.Circle(centrecircle,r,fill=False,
                                            edgecolor=circle_settings['edgecolor'],
                                            linestyle=circle_settings['linestyle'],
                                            linewidth=circle_settings['linewidth'],
                                            alpha=0.5)
                        ax[i,j].add_patch(circle)

        hdul.close()
    fig.subplots_adjust(hspace=0,wspace=0,left=0,right=1, bottom=0, top=1)
    if returnfig:
        return fig
    else:
        fig.savefig('plot_comp_facevsedge.png',format='png',dpi=300)


def plot_single_var_projection(proj_FITS,field,logscale=True,scalebar=(3,'kpc'),redshift=True,colorbar=True,
                                colormap=None, type_scale='',centers=[],radii=[],names=[],filter_name='none',
                                smooth=False):
    """This function uses the projection information in a FITS file following the 
        OZY format and plots it following the OZY standards."""
    
    # Make required imports
    from ozy.utils import invert_tick_colours
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    from matplotlib.colors import LogNorm,SymLogNorm
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import seaborn as sns
    from astropy.convolution import Gaussian2DKernel, convolve
    import scipy as sp   
    sns.set(style="dark")
    plt.rcParams["axes.axisbelow"] = False
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    # hfont = {'fontname':'Helvetica'}
    # matplotlib.rc('text', usetex = True)
    # matplotlib.rc('font', **{'family' : "serif"})
    # params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
    # matplotlib.rcParams.update(params)

    stellar = False
    # First,check that FITS file actually exists
    if not os.path.exists(proj_FITS):
        raise ImportError('File not found. Please check!')
    
    # Load FITS file
    hdul = fits.open(proj_FITS)
    hdul_fields = [h.header['btype'] for h in hdul]

    # Add filter identification to the field
    hdul_filtername = ['COND_0' in h.header for h in hdul]
    if any(hdul_filtername):
            field = field + '/' + filter_name
    else:
        print('WARNING: This FITS file does not have filter information, so anything you added will be ignored!')

    # Check that the required fields for plotting are in this FITS
    if field not in hdul_fields:
        print('Field %s is not between the ones you provided... Check!'%field)
        exit

    if type_scale != '':
        type_scale = '_'+type_scale

    # Since everything is fine, we begin plotting…
    figsize = plt.figaspect(float(7) / float(7))
    fig = plt.figure(figsize=figsize, facecolor='k', edgecolor='k')
    fig, ax = plt.subplots(1, 1, figsize=(7,7), dpi=200, facecolor='k', edgecolor='k')

    width_x =  hdul[0].header['CDELT1']*hdul[0].header['NAXIS1']
    width_y =  hdul[0].header['CDELT2']*hdul[0].header['NAXIS2']
    ex = [-0.5*width_x,0.5*width_x,-0.5*width_y,0.5*width_y]
    ax.set_xlim([-0.5*width_x,0.5*width_y])
    ax.set_ylim([-0.5*width_x,0.5*width_y])
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ax.axis('off')
    h = [k for k in range(0,len(hdul)) if hdul[k].header['btype']==field][0]
    if field.split('/')[0] == 'star' or field.split('/')[0] == 'dm':
        plotting_def = get_plotting_def(field.split('/')[0]+'_'+field.split('/')[1])
        stellar = True
        full_varname = field.split('/')[0]+'_'+field.split('/')[1]
    else:
        plotting_def = get_plotting_def(field.split('/')[1])
        full_varname = field.split('/')[1]
    if colormap == None:
        colormap = plotting_def['cmap']
    if smooth:
        # We smooth with a Gaussian kernel with x_stddev=1 (and y_stddev=1)
        # It is a 9x9 array
        kernel = Gaussian2DKernel(x_stddev=0.8)
        # astropy's convolution replaces the NaN pixels with a kernel-weighted
        # interpolation from their neighbors
        cImage = convolve(hdul[h].data.T, kernel)
        # sigma=3
        # cImage = sp.ndimage.filters.gaussian_filter(hdul[h].data.T, sigma, mode='constant')
    else:
        cImage = hdul[h].data.T
    print(field,np.nanmin(cImage),np.nanmax(cImage))
    if logscale and not plotting_def['symlog']:
        plot = ax.imshow(np.log10(cImage), cmap=colormap,
                        origin='upper',vmin=np.log10(plotting_def['vmin'+type_scale]),
                        vmax=np.log10(plotting_def['vmax'+type_scale]),extent=ex,
                        interpolation='nearest')
    elif logscale and not plotting_def['symlog']:
        if (filter_name != 'outflow' and filter_name != 'inflow'):
            plot = ax.imshow(cImage, cmap=colormap,
                            origin='upper',norm=SymLogNorm(linthresh=0.1, linscale=1,
                            vmin=plotting_def['vmin'+type_scale], vmax=plotting_def['vmax'+type_scale]),
                            extent=ex,
                            interpolation='nearest')
        elif field.split('/')[1] == 'v_sphere_r':
            if filter_name == 'inflow' and smooth:                                                                                                                                                                                          
                cImage = sp.ndimage.filters.gaussian_filter(-hdul[h].data.T, sigma, mode='constant')
            plot = ax.imshow(np.log10(cImage), cmap=plotting_def['cmap_'+filter_name],
                            origin='upper',vmin=np.log10(plotting_def['vmin_'+filter_name]),
                            vmax=np.log10(plotting_def['vmax_'+filter_name]),extent=ex,
                            interpolation='nearest')
        else:
            plot = ax.imshow(cImage, cmap=colormap,
                            origin='upper',norm=SymLogNorm(linthresh=0.1, linscale=1,
                            vmin=plotting_def['vmin'+type_scale], vmax=plotting_def['vmax'+type_scale]),
                            extent=ex,
                            interpolation='nearest')
    else:
        plot = ax.imshow(hdul[h].data.T, cmap=colormap,
                        origin='upper',extent=ex,interpolation='nearest',
                        vmin=plotting_def['vmin'+type_scale],vmax=plotting_def['vmax'+type_scale])
    if colorbar:
        cbaxes = inset_axes(ax, width="80%", height="5%", loc='lower center')
        cbar = fig.colorbar(plot, cax=cbaxes, orientation='horizontal')
        if logscale and field.split('/')[1] != 'v_sphere_r':
            cbar.set_label(plotting_def['label_log'],color=plotting_def['text_over'],fontsize=20,labelpad=-60, y=0.85,weight='bold')
        else:
            cbar.set_label(plotting_def['label'],color=plotting_def['text_over'],fontsize=20,labelpad=-10, y=1.25)
        cbar.ax.xaxis.label.set_font_properties(matplotlib.font_manager.FontProperties(weight='bold',size=15))
        cbar.ax.tick_params(axis='x', pad=-16, labelsize=13,labelcolor=plotting_def['text_over'])
        cbar.ax.tick_params(length=0,width=0)
        invert_tick_colours(cbar.ax,full_varname,type_scale)

    if redshift:
        ax.text(0.03, 0.95, r'$z = {z:.2f}$'.format(z=hdul[h].header['redshift']), # Redshift
                            verticalalignment='bottom', horizontalalignment='left',
                            transform=ax.transAxes,
                            color=plotting_def['text_over'], fontsize=20,fontweight='bold')

    fontprops = fm.FontProperties(size=20,weight='bold')
    if scalebar != None:
        sl = int(unyt_quantity(scalebar[0],scalebar[1]).in_units(hdul[0].header['CUNIT1']).d)
        slb = AnchoredSizeBar(ax.transData,
                                    sl, '%s %s'%(str(int(sl)),hdul[0].header['CUNIT1']), 
                                    'upper right', 
                                    pad=0.1,
                                    color=plotting_def['text_over'],
                                    frameon=False,
                                    size_vertical=0.1,
                                    fontproperties=fontprops)
        ax.add_artist(slb)

    if len(centers) != 0 and len(radii) != 0:
        for c in range(0, len(centers)):
            centrecircle = (-centers[c][2]*1000,-centers[c][1]*1000)
            r = radii[c] * 1000
            circle = plt.Circle(centrecircle,r,fill=False,edgecolor='w',linestyle='--')
            ax.add_patch(circle)
            ax.text(centrecircle[0]+1.1*r, centrecircle[1]+1.1*r, names[c], # Name of object
                            verticalalignment='bottom', horizontalalignment='left',
                            color=plotting_def['text_over'], fontsize=10,fontweight='bold')
        # centrecircle = (-centers[0][2]*1000,-centers[0][1]*1000)
        # r = radii[0] * 1000
        # circle = plt.Circle(centrecircle,r,fill=False,edgecolor='w',linestyle='--')
        # ax.add_patch(circle)
        # ax.text(centrecircle[0]+1.1*r, centrecircle[1]+1.1*r, names[0], # Name of object
        #                 verticalalignment='bottom', horizontalalignment='left',
        #                 color=plotting_def['text_over'], fontsize=10,fontweight='bold')
        # ax.scatter(-centers[1:][2]*1000,-centers[1:][1]*1000)#,s=0.1,alpha=0.4,facecolor='r')

    fig.subplots_adjust(hspace=0,wspace=0,left=0,right=1, bottom=0, top=1)

    if filter_name != 'none':
        fig.savefig(proj_FITS.split('.fits')[0]+'_'+field.split('/')[1]+'_'+filter_name+'.png',format='png',dpi=330)
    else:
        fig.savefig(proj_FITS.split('.fits')[0]+'_'+field.split('/')[1]+'.png',format='png',dpi=330)

def plot_lupton_rgb_projection(proj_FITS,fields,stars=False,scalebar=(3,'kpc'),redshift=True, type_scale='galaxy'):
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
            plotting_def = get_plotting_def(fields[i].split('/')[0]+'_'+fields[i].split('/')[1])
        else:
            plotting_def = get_plotting_def(fields[i].split('/')[1])

        data = hdul[h].data.T
        print(fields[i],data.min(),data.max())
        if type_scale == 'galaxy':
            data[data < plotting_def['vmin'+type_scale]] = plotting_def['vmin'+type_scale]
            data[data > plotting_def['vmax'+type_scale]] = plotting_def['vmax'+type_scale]
            data = (np.log10(data)-np.log10(plotting_def['vmin'+type_scale]*0.99))/(np.log10(plotting_def['vmax'+type_scale]) - np.log10(plotting_def['vmin'+type_scale]*0.99))
        else:
            data[data < plotting_def['vmin']] = plotting_def['vmin']
            data[data > plotting_def['vmax']] = plotting_def['vmax']
            data = (np.log10(data)-np.log10(plotting_def['vmin']*0.99))/(np.log10(plotting_def['vmax']) - np.log10(plotting_def['vmin']*0.99))
        images.append(data)
        custom_lines.append(mpatches.Patch(color=rgb_colours[i],label=plotting_def['label']))
    
    rgb_default = make_lupton_rgb(images[0], images[1], images[2],Q=10, stretch=0.5)

    ax.imshow(rgb_default,origin='upper',extent=ex,interpolation='nearest')

    if stars:
        h = [k for k in range(0,len(hdul)) if hdul[k].header['btype']=='star/mass'][0]
        stars = np.log10(hdul[h].data.T)
        plotting_def = get_plotting_def('star_mass')
        ax.imshow(stars,origin='upper',cmap=plotting_def['cmap'],
                    extent=ex,interpolation='nearest',
                    vmin=np.log10(plotting_def['vmin']),
                    vmax=np.log10(plotting_def['vmax']),alpha=0.4)

    if redshift:
        ax.text(0.05, 0.90, r'$z = {z:.2f}$'.format(z=hdul[h].header['redshift']), # Redshift
                    verticalalignment='bottom', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='white', fontsize=16,fontweight='bold')

    fontprops = fm.FontProperties(size=16,weight='bold')
    if scalebar != None:
        sl = unyt_quantity(scalebar[0],scalebar[1]).in_units(hdul[0].header['CUNIT1']).d
        scalebar = AnchoredSizeBar(ax.transData,
                                    sl, '%s %s'%(str(sl),hdul[0].header['CUNIT1']),
                                    'upper right', 
                                    pad=0.1,
                                    color=plotting_def['text_over'],
                                    frameon=False,
                                    size_vertical=0.1,
                                    fontproperties=fontprops)
        ax.add_artist(scalebar)
    ax.legend(handles=custom_lines,loc='lower right',frameon=False,fontsize=12,labelcolor='white')
    fig.subplots_adjust(hspace=0,wspace=0,left=0,right=1, bottom=0, top=1)
    fig.savefig(proj_FITS.split('.')[0]+'_rgb.png',format='png',dpi=300)


def do_healpix_projection(group,vars,weight=['gas/density','star/age'],nside=32,pov='edgeon',
                          r=(1.0,'rvir'),dr=(1./150.,'rvir'),rmin=None,lmax=100,lmin=1,remove_subs=False,
                          filter_conds=['none'],filter_name=['none'],use_neigh=False,myaxis=np.array([1.,0.,0.])):
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

    # Setup camera details for the requested POV (Point of View)
    centre = vectors.vector()
    centre.x, centre.y, centre.z = group.position[0], group.position[1], group.position[2]
    bulk = vectors.vector()
    velocity = group.velocity.in_units('code_velocity')
    bulk.x, bulk.y, bulk.z = velocity.d[0], velocity.d[1], velocity.d[2]
    region_axis = vectors.vector()
    norm_L = group.angular_mom['total'].d/np.linalg.norm(group.angular_mom['total'].d)
    region_axis.x,region_axis.y,region_axis.z = norm_L[0], norm_L[1], norm_L[2]
    
    if pov == 'edgeon':
        los = cartesian_basis['x'] - np.dot(cartesian_basis['x'],norm_L)*norm_L
        los /= np.linalg.norm(los)
        axis = vectors.vector()
        axis.x,axis.y,axis.z = los[0], los[1], los[2]
        up_vector = vectors.vector()
        up_vector.x,up_vector.y,up_vector.z = norm_L[0], norm_L[1], norm_L[2]
        if r[1] == 'rvir':
            rmax = r[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius']
        else:
            rmax = group.obj.quantity(r[0],r[1]).in_units('code_length')
        if rmin == None:
            if dr[1] == 'rvir':
                dr = dr[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius']
            else:
                dr = group.obj.quantity(dr[0],dr[1]).in_units('code_length')
            distance = rmax - 0.5*dr
            far_cut_depth = rmax + 0.5*dr
        else:
            if rmin[1] == 'rvir':
                rmin = rmin[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius']
            else:
                rmin = group.obj.quantity(rmin[0],rmin[1]).in_units('code_length')
            distance = rmin
            far_cut_depth = rmax
        region_size = np.array([2.0*far_cut_depth,2.0*far_cut_depth],order='F',dtype=np.float64)

        # Update projection details with the ones used for the region
        proj.up_vector = norm_L
        proj.centre = group.position[:]
        proj.nside = nside
    elif pov == 'custom':
        norm_L = myaxis/np.linalg.norm(myaxis)
        los = cartesian_basis['x'] - np.dot(cartesian_basis['x'],norm_L)*norm_L
        los /= np.linalg.norm(los)
        axis = vectors.vector()
        axis.x,axis.y,axis.z = los[0], los[1], los[2]
        up_vector = vectors.vector()
        up_vector.x,up_vector.y,up_vector.z = norm_L[0], norm_L[1], norm_L[2]
        if r[1] == 'rvir':
            rmax = r[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius']
        else:
            rmax = group.obj.quantity(r[0],r[1]).in_units('code_length')
        if rmin == None:
            if dr[1] == 'rvir':
                dr = dr[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius']
            else:
                dr = group.obj.quantity(dr[0],dr[1]).in_units('code_length')
            distance = rmax - 0.5*dr
            far_cut_depth = rmax + 0.5*dr
        else:
            if rmin[1] == 'rvir':
                rmin = rmin[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius']
            else:
                rmin = group.obj.quantity(rmin[0],rmin[1]).in_units('code_length')
            distance = rmin
            far_cut_depth = rmax
        region_size = np.array([2.0*far_cut_depth,2.0*far_cut_depth],order='F',dtype=np.float64)

        # Update projection details with the ones used for the region
        proj.up_vector = norm_L
        proj.centre = group.position[:]
        proj.nside = nside
    else:
        print('This POV for a HEALPix projection is not supported...')
        print('Please check!')
        exit
        
    # Now create filters if any conditions have been given...
    nfilter = len(filter_conds)
    for i in range(0,nfilter):
        f = init_filter(filter_conds[i],filter_name[i],group)
        proj._get_python_filter(f)
        
    # Do it for hydro first
    filts = []
    for i in range(0,nfilter):
        if isinstance(filter_conds[i],list):
            cond_var = filter_conds[i][0].split('/')[0]
        else:
            cond_var = filter_conds[i].split('/')[0]
        if cond_var in common_variables or cond_var in grid_variables:
            f = init_filter(filter_conds[i],filter_name[i],group)
        else:
            # When a filter asks for a variable not existent in the common_variables
            # or the grid_variables dictionaries just ignore it and set it to blank
            f = init_filter('none','none',group)
        filts.append(f)
        
    # Construct substructure regions if we want them out of the projection
     # Construct substructure regions if we want them out of the projection
    remove_all = False
    if remove_subs == 'all':
        remove_all = True
        remove_subs = True
    if remove_all:
        subs = structure_regions(group, add_substructure=False, add_neighbours=False,
                                    add_intersections=True, position=proj.centre,
                                    radius=far_cut_depth)
        nsubs = len(subs)
    elif remove_subs:
        subs = structure_regions(group, add_substructure=True, add_neighbours=False,
                                    tidal_method='BT87_simple')
        nsubs = len(subs)
    else:
        nsubs = 0
        
    cam = obs_instruments.init_camera(centre,axis,up_vector,region_size,region_axis,bulk,distance,
                                      far_cut_depth,1024,nfilter,nsubs)
    
    # Now give filters to camera Fortran type
    for i in range(0,nfilter):
        cam.filters[i] = filts[i]
        
    if remove_subs and nsubs>0:
        for i in range(0,nsubs):
            cam.subs[i] = subs[i]

    
    # Create projection_handler Fortran derived type for the results of the hydro data projection
    hydro_handler = maps.projection_handler()
    hydro_handler.type = pov
    hydro_handler.nvars = len(proj.vars['gas'])
    hydro_handler.weightvar = proj.weight[0]
    maps.allocate_projection_handler(hydro_handler)
    for i in range(0, len(proj.vars['gas'])):
        hydro_handler.varnames.T.view('S128')[i] = proj.vars['gas'][i].ljust(128)
    
    # COMPUTE HYDRO PROJECTION
    if lmax != 0 or lmin != 0:
        maps.healpix_hydro(group.obj.simulation.fullpath,cam,use_neigh,hydro_handler,nside,int(lmax),int(lmin))
    else:
        maps.healpix_hydro(group.obj.simulation.fullpath,cam,use_neigh,hydro_handler,nside)
    # TODO: Weird issue when the direct map array is given.
    data = np.copy(hydro_handler.map)
    proj.data_maps.append(data)
    # TODO: Add particle projections

    return proj
    
def load_single_galaxy_healpix(proj_FITS,field,filter_name='none'):
    """This function uses the HEALPix projection information in a FITS file following
        the OZY format and returns the map as an array."""
    
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
    plt.rcParams["axes.axisbelow"] = False
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    # hfont = {'fontname':'Helvetica'}
    # matplotlib.rc('text', usetex = True)
    # matplotlib.rc('font', **{'family' : "serif"})
    # params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
    # matplotlib.rcParams.update(params)

    # First, check that the FITS file actually exists
    if not os.path.exists(proj_FITS):
        raise ImportError('File not found. Please check!')
    
    # Load FITS file
    hdul = fits.open(proj_FITS)
    hdul_fields = hdul[1].data.columns.names
        
    # Add filter identification to the field
    hdul_filtername = ['COND_0' in h.header for h in hdul]
    if any(hdul_filtername):
        field = field + '/' + filter_name
    else:
        print('WARNING: This FITS file does not have filter information, so anything you added will be ignored!')
        
    # Check that the required field is in this FITS
    if field not in hdul_fields:
        print('The field %s is not included in this file. Ignoring...'%field)
        exit
       
    # Make adaptations for HEALPix and Astropy
    target_wcs = WCS(target_header)

    # Since everything is fine, we begin plotting…
    nside = hdul[1].header["NSIDE"]

    data = hdul[1].data[field]

    return data,nside,target_wcs

def plot_single_galaxy_healpix(proj_FITS,fields,logscale=True,redshift=False,fig_format='png',
                               type_scale='',filter_name='none'):
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
    from ozy.utils import invert_tick_colours
    sns.set_theme(style="dark")
    plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "Computer Modern Roman",
    })

    # First, check that the FITS file actually exists
    if not os.path.exists(proj_FITS):
        raise ImportError('File not found. Please check!')
    
    # Add filter identification to the field
    hdul = fits.open(proj_FITS)
    hdul_filtername = ['COND_0' in h.header for h in hdul]
    if any(hdul_filtername):
        for f in range(0,len(fields)):
            fields[f] = fields[f] + '/' + filter_name
    else:
        print('WARNING: This FITS file does not have filter information, so anything you added will be ignored!')
    
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
    # Since everything is fine, we begin plotting…
    nside = hdul.header["NSIDE"]
    ncolumns = int(len(fields)/2)
    figsize = plt.figaspect(float(6 * 2) / float(9 * ncolumns))
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
            data[invalid_index] = data[invalid_index-1]
            array, footprint = reproject_from_healpix((data, 'C'),
                                                        target_header,nested=False)
            print(fields[ivar],np.nanmin(array),np.nanmax(array))
            
            if fields[ivar].split('/')[0] == 'star' or fields[ivar].split('/')[0] == 'dm':
                plotting_def = get_plotting_def(fields[ivar].split('/')[0]+'_'+fields[ivar].split('/')[1])
                stellar = True
            else:
                plotting_def = get_plotting_def(fields[ivar].split('/')[1])
            if logscale:
                if not plotting_def['symlog']:
                    plot = ax.imshow(np.log10(array), cmap=plotting_def['cmap'],
                                    vmin=np.log10(plotting_def['vmin'+type_scale]),
                                    vmax=np.log10(plotting_def['vmax'+type_scale]))
                else:
                    plot = ax.imshow(array, cmap=plotting_def['cmap'],
                                    norm=SymLogNorm(linthresh=plotting_def['linthresh'], 
                                                                   linscale=plotting_def['linscale'],
                                                                   vmin=plotting_def['vmin'+type_scale], 
                                                                   vmax=plotting_def['vmax'+type_scale]))
            else:
                plot = ax.imshow(array, cmap=plotting_def['cmap'],
                                vmin=plotting_def['vmin'+type_scale],vmax=plotting_def['vmax'+type_scale])
            ax.coords.grid(color=plotting_def['text_over'],linewidth=0.8,alpha=0.4,linestyle=':')
            # ax.coords['ra'].set_ticklabel(color=plotting_def['text_over'])
            ax.coords['ra'].set_ticklabel_visible(False)
            ax.coords['ra'].set_ticks_visible(False)
            ax.coords['dec'].set_ticks_position('lb')
            ax.coords['dec'].set_ticklabel_position('lb')
            fontprops = fm.FontProperties(size=14)
            if i == 0:
                cbaxes = inset_axes(ax, width="80%", height="10%", loc='upper center',
                                    bbox_to_anchor=(0.0, 0.2, 1, 1.),
                                    bbox_transform=ax.transAxes)
                axcb = fig.colorbar(plot, cax = cbaxes, orientation='horizontal',pad=0.3)
                if logscale and not plotting_def['symlog']:
                    axcb.set_label(plotting_def['label_log'],color='black',fontsize=14,labelpad=-32)
                else:
                    axcb.set_label(plotting_def['label'],color='black',fontsize=14,labelpad=-32)
                axcb.ax.tick_params(axis='x', pad=-10, labelsize=12,labelcolor='white')
                axcb.outline.set_linewidth(1)
                axcb.ax.set_axisbelow(False)
                invert_tick_colours(axcb.ax,fields[ivar].split('/')[1],'_cgm')
                axcb.ax.tick_params(length=0,width=0)
            else:
                cbaxes = inset_axes(ax, width="80%", height="10%", loc='lower center',
                                    bbox_to_anchor=(0.0, -0.2, 1., 1.),
                                    bbox_transform=ax.transAxes)
                axcb = fig.colorbar(plot, cax = cbaxes, orientation='horizontal',pad=0.3)
                if logscale and not plotting_def['symlog']:
                    axcb.set_label(plotting_def['label_log'],color='black',fontsize=14)
                else:
                    axcb.set_label(plotting_def['label'],color='black',fontsize=14)
                axcb.ax.tick_params(axis='x', pad =-10, labelsize=12,labelcolor='white')
                axcb.outline.set_linewidth(1)
                axcb.ax.set_axisbelow(False)
                invert_tick_colours(axcb.ax,fields[ivar].split('/')[1],'_cgm')
                axcb.ax.tick_params(length=0,width=0)
    if redshift:
        fig.text(0.5, 0.5, 'z = '+str(round(hdul.header['REDSHIFT'], 2)),
                fontsize=18,va='center',
                color=plotting_def['text_over'])
    fig.subplots_adjust(hspace=0, wspace=0,top = 0.90,bottom = 0.12,left = 0.02,right = 0.98)
    if filter_name != 'none':
        fig.savefig(proj_FITS.split('.')[0]+'_'+filter_name+f'.{fig_format}',format=fig_format,dpi=300)
    else:
        fig.savefig(proj_FITS.split('.')[0]+f'.{fig_format}',format=fig_format,dpi=300)

def structure_overplotting(group, add_substructure=True, add_neighbours=False,
                            tidal_method='BT87_simple',rmax=(1e10,'kpc')):
    """
    This routine obtains the plotting details to speed up the construction of
    projections with substructure and/or neighbours information.
    """
    from ozy.variables_settings import circle_dictionary
    from ozy.utils import tidal_radius
    # Setup arrays with information
    centers = []
    radii = []
    plot_dicts = []

    # Get the central object information
    if group.type == 'halo':
        myhalo = group
        centers.append(group.position)
        radii.append(myhalo.virial_quantities['radius'])
        plot_dicts.append(circle_dictionary['rvir_halo'])
    elif group.type == 'galaxy':
        myhalo = group.halo
        centers.append(group.position)
        radii.append(0.2*myhalo.virial_quantities['radius'])
        plot_dicts.append(circle_dictionary['rvir_galaxy'])

    if isinstance(rmax, tuple):
        rmax = group.obj.quantity(rmax[0],rmax[1])
    # If asked for substructure, obtain the substructure of the host halo
    if add_substructure:
        subs = myhalo.substructure_list
        for s in subs:
            # Get halo galaxies
            sub_gals = s.galaxies
            position = s.position
            mysub = s
            if len(sub_gals) != 0:
                for sg in sub_gals:
                    if sg.central:
                        position = sg.position
                        mysub = sg
            distance = group.position - position
            d = group.obj.quantity(np.linalg.norm(distance.to('kpc').d),'kpc')
            if s.npart >= 1000 and d.to('kpc')<=rmax.to('kpc'):
                centers.append(mysub.position)
                radii.append(s.radius[tidal_method])
                plot_dicts.append(circle_dictionary['tidal_'+tidal_method])
    
    # If asked for neighbours (so inside the virial radius) obtain them
    if add_neighbours:
        # TODO: This needs to be updated with the actual class procedure
        # myhalo.get_neighbours_in(rvir)
        pass
    return centers,radii,plot_dicts