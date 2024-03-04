from matplotlib.pyplot import plot
import numpy as np
from collections import deque,Counter
from bisect import insort, bisect_left
from itertools import islice
import sys
import os
import subprocess
from ozy.dict_variables import get_code_units
import matplotlib.text as mtext
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt

class RotationAwareAnnotation(mtext.Annotation):
    def __init__(self, s, xy, p, pa=None, ax=None, **kwargs):
        self.ax = ax or plt.gca()
        self.p = p
        if not pa:
            self.pa = xy
        kwargs.update(rotation_mode=kwargs.get("rotation_mode", "anchor"))
        mtext.Annotation.__init__(self, s, xy, **kwargs)
        self.set_transform(mtransforms.IdentityTransform())
        if 'clip_on' in kwargs:
            self.set_clip_path(self.ax.patch)
        self.ax._add_text(self)

    def calc_angle(self):
        p = self.ax.transData.transform_point(self.p)
        pa = self.ax.transData.transform_point(self.pa)
        ang = np.arctan2(p[1]-pa[1], p[0]-pa[0])
        return np.rad2deg(ang)

    def _get_rotation(self):
        return self.calc_angle()

    def _set_rotation(self, rotation):
        pass

    _rotation = property(_get_rotation, _set_rotation)

def most_contrast_rgba(rgba):
    """
    Returns the most contrasting RGBA color for a given RGBA color.
    """
    # Extract the RGBA components
    red, green, blue, alpha = rgba

    # Calculate the luminance of the color
    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue

    # Calculate the opposite color
    opposite_red = 1 - red
    opposite_green = 1 - green
    opposite_blue = 1 - blue

    # Calculate the opposite color's luminance
    opposite_luminance = 0.2126 * opposite_red + 0.7152 * opposite_green + 0.0722 * opposite_blue

    return (opposite_red, opposite_green, opposite_blue, alpha)
    # # If the luminance of the opposite color is greater, return the opposite color
    # if opposite_luminance > luminance:
    #     return (opposite_red, opposite_green, opposite_blue, alpha)

    # # Otherwise, return black or white depending on the luminance of the original color
    # if luminance < 0.5:
    #     return (0, 0, 0, alpha) # Black
    # else:
    #     return (1, 1, 1, alpha) # White


def invert_tick_colours(ax,var,type_scale):
    from plot_settings import plotting_dictionary, symlog_variables
    from matplotlib.colors import LogNorm,SymLogNorm
    from matplotlib import colormaps

    fig = plt.gcf()
    plotting_def = plotting_dictionary[var]
    cmap = colormaps.get_cmap(plotting_def['cmap'])
    ticks_pos = ax.get_xticks()
    ticks_labels = ax.get_xticklabels()
    if var not in symlog_variables:
        norm = LogNorm(vmin=plotting_def['vmin'+type_scale],
                         vmax=plotting_def['vmax'+type_scale],
                         clip=True)
        for tp,tl in zip(ticks_pos,ticks_labels):
            rgba = cmap(norm(10**tp))
            new_rgba = most_contrast_rgba(rgba)
            tl.set_color(new_rgba)
        fig.canvas.draw()
    else:
        norm = SymLogNorm(vmin=plotting_def['vmin'+type_scale],
                         vmax=plotting_def['vmax'+type_scale],
                         linthresh=plotting_def['linthresh'],
                         linscale=plotting_def['linscale'],
                         clip=True)
        for tp,tl in zip(ticks_pos,ticks_labels):
            rgba = cmap(norm(tp))
            new_rgba = most_contrast_rgba(rgba)
            tl.set_color(new_rgba)
        fig.canvas.draw()
        

def as_si(x, ndp):
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))
    
def read_infofile(infopath):
    info = {}
    info['aexp'] = 0.0
    info['redshift'] = 0.0

    # Open info file and read the units conversion factors from it
    with open(infopath,'r') as infofile:
        cc = 1
        while (cc <= 20):
            newline = infofile.readline()
            linehead = newline[0:12].strip()
            if (linehead == "ncpu")  : info['ncpu']  = int(newline[13:].strip())
            if (linehead == "ndim")  : info['ndim'] = int(newline[13:].strip())
            if (linehead == "levelmin")  : info['levelmin']  = int(newline[13:].strip())
            if (linehead == "levelmax")  : info['levelmax'] = int(newline[13:].strip())
            if (linehead == "ngridmax")  : info['ngridmax']  = int(newline[13:].strip())
            if (linehead == "nstep_coarse")  : info['nstep_coarse'] = int(newline[13:].strip())

            if (linehead == "aexp")  : info['aexp']  = float(newline[13:].strip())
            if (linehead == "time")  : info['time'] = float(newline[13:].strip())
            if (linehead == "H0"): info['H0']  = float(newline[13:].strip())
            if (linehead == "omega_m"): info['omega_m']  = float(newline[13:].strip())
            if (linehead == "omega_l"): info['omega_l']  = float(newline[13:].strip())
            if (linehead == "omega_k"): info['omega_k']  = float(newline[13:].strip())
            if (linehead == "omega_b"): info['omega_b']  = float(newline[13:].strip())
            if (linehead == "unit_l"): info['unit_l']  = float(newline[13:].strip())
            if (linehead == "unit_d"): info['unit_d']  = float(newline[13:].strip())
            if (linehead == "unit_t"): info['unit_t']  = float(newline[13:].strip())

            cc = cc + 1
        
        if info['aexp'] < 1.0 and info['time'] <= 0.0:
            info['redshift'] = 1./info['aexp'] - 1.
    return info

def closest_snap_z(simfolder,z,return_index=False):
    """
    Using the IDtoZetas Fortran script, it gets the snapshot
    closest in redshift to the wanted value.
    """
    presentpath = os.getcwd()
    os.chdir(simfolder)
    result = subprocess.run('IDtoZetas.out -ask '+str(z), shell=True,stdout=subprocess.PIPE)
    os.chdir(presentpath)
    indexout = int(result.stdout.decode('utf-8').split('is')[1].split('(z')[0])

    ozyfile = 'ozy_%05d.hdf5' % (indexout)

    if return_index:
        return ozyfile, indexout
    else:
        return ozyfile
    
def get_tdyn(galaxy):
        """
        Computes the dynamical time-scale tdyn as
        the time required for a test particle to complete
        one full orbit at 0.2 Rvir.

        tdyn = 2pi*sqrt(R^3/(GM))
        where we assume M = Mgas+Mstars*Mdm is measured 
        within 0.2 Rvir
        """
        from unyt import G
        

        Mtot = galaxy.mass['dm'] + galaxy.mass['baryon']
        r = 0.2*galaxy.obj.halos[galaxy.parent_halo_index].virial_quantities['radius']
        tdyn = 2*np.pi*np.sqrt(r**3/(G*Mtot))
        return tdyn

def find_neigh_snaps(simfolder,orig_snap,trange,minsnaps=3,returnweight=False):
    """
    This function searches for the closest snapshots
    to a given one within a particular time frame.

    It also gives the option of returning weights
    corresponding to the contribution of each one to the time considered.
    TODO: Use 3 snaps as a minimum
    """
    import glob
    from astropy.cosmology import FlatLambdaCDM
    import ozy

    # Firstly, list all snapshots in the directory
    presentpath = os.getcwd()
    os.chdir(simfolder)
    snapshots = glob.glob('output_0*')
    snapshots.sort(key=lambda x: int(x[-5:]))
    iorig = snapshots.index(orig_snap)

    # Get time of original snapshot
    ozy_orig = 'ozy_%05d.hdf5'%(int(snapshots[iorig][-5:]))
    sim = ozy.load('Groups/'+ozy_orig)
    cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                                    Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)
    t_orig = cosmo.age(sim.simulation.redshift).value

    neigh_snaps = []
    weights = []
    times = []
    # Find the snapshots just below the original one
    for i in range(iorig-1,0,-1):
        ozy_name = 'ozy_%05d.hdf5'%(int(snapshots[i][-5:]))
        sim = ozy.load('Groups/'+ozy_name)
        cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                                        Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)
        thubble = cosmo.age(sim.simulation.redshift).value
        try:
            next_ozy_name = 'ozy_%05d.hdf5'%(int(snapshots[i-1][-5:]))
            next_sim = ozy.load('Groups/'+next_ozy_name)
            next_cosmo = FlatLambdaCDM(H0=next_sim.simulation.hubble_constant, Om0=next_sim.simulation.omega_matter, 
                                            Ob0=next_sim.simulation.omega_baryon,Tcmb0=2.73)
            t_next = cosmo.age(next_sim.simulation.redshift).value
        except:
            t_next = t_orig - 0.5*trange
        prev_ozy_name = 'ozy_%05d.hdf5'%(int(snapshots[i+1][-5:]))
        prev_sim = ozy.load('Groups/'+prev_ozy_name)
        prev_cosmo = FlatLambdaCDM(H0=prev_sim.simulation.hubble_constant, Om0=prev_sim.simulation.omega_matter, 
                                        Ob0=prev_sim.simulation.omega_baryon,Tcmb0=2.73)
        t_prev = cosmo.age(prev_sim.simulation.redshift).value
        if t_orig - thubble <= 0.5*trange and t_orig > thubble:
            neigh_snaps.append(ozy_name)
            tup = 0.5*(t_prev - thubble)
            tdown = 0.5*(thubble - t_next)
            weights.append((tup+tdown)/abs(t_orig-thubble))
            times.append(thubble)
        elif t_orig - thubble > 0.5*trange:
            if len(neigh_snaps) == 0:
                # In the case that we need to extend a bit further
                neigh_snaps.append(ozy_name)
                tup = 0.5*(t_prev - thubble)
                weights.append(tup)
                times.append(thubble)
                trange = 2*tup
            break
    
    # Add original snapshot
    neigh_snaps.append(ozy_orig)
    next_ozy_name = 'ozy_%05d.hdf5'%(int(snapshots[iorig-1][-5:]))
    next_sim = ozy.load('Groups/'+next_ozy_name)
    next_cosmo = FlatLambdaCDM(H0=next_sim.simulation.hubble_constant, Om0=next_sim.simulation.omega_matter, 
                                    Ob0=next_sim.simulation.omega_baryon,Tcmb0=2.73)
    t_next = cosmo.age(next_sim.simulation.redshift).value
    prev_ozy_name = 'ozy_%05d.hdf5'%(int(snapshots[iorig+1][-5:]))
    prev_sim = ozy.load('Groups/'+prev_ozy_name)
    prev_cosmo = FlatLambdaCDM(H0=prev_sim.simulation.hubble_constant, Om0=prev_sim.simulation.omega_matter, 
                                    Ob0=prev_sim.simulation.omega_baryon,Tcmb0=2.73)
    t_prev = cosmo.age(prev_sim.simulation.redshift).value
    # And compute the weight of the original/middle snapshot
    tup = 0.5*(t_prev - t_orig)
    tdown = 0.5*(t_orig - t_next)
    weights.append(tup+tdown)
    times.append(t_orig)
    indexorig = len(weights)
    # And do the same for just above the original one
    for i in range(iorig+1,len(snapshots),1):
        ozy_name = 'ozy_%05d.hdf5'%(int(snapshots[i][-5:]))
        sim = ozy.load('Groups/'+ozy_name)
        cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                                        Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)
        thubble = cosmo.age(sim.simulation.redshift).value
        try:
            prev_ozy_name = 'ozy_%05d.hdf5'%(int(snapshots[i+1][-5:]))
            prev_sim = ozy.load('Groups/'+prev_ozy_name)
            prev_cosmo = FlatLambdaCDM(H0=prev_sim.simulation.hubble_constant, Om0=prev_sim.simulation.omega_matter, 
                                            Ob0=prev_sim.simulation.omega_baryon,Tcmb0=2.73)
            t_prev = cosmo.age(prev_sim.simulation.redshift).value
        except:
            t_prev = t_orig + 0.5*trange
        next_ozy_name = 'ozy_%05d.hdf5'%(int(snapshots[i-1][-5:]))
        next_sim = ozy.load('Groups/'+next_ozy_name)
        next_cosmo = FlatLambdaCDM(H0=next_sim.simulation.hubble_constant, Om0=next_sim.simulation.omega_matter, 
                                        Ob0=next_sim.simulation.omega_baryon,Tcmb0=2.73)
        t_next = cosmo.age(next_sim.simulation.redshift).value
        if thubble - t_orig <= 0.5*trange and t_orig < thubble:
            neigh_snaps.append(ozy_name)
            tdown = 0.5*(thubble - t_next)
            tup = 0.5*(t_prev - thubble)
            weights.append((tup+tdown)/abs(t_orig-thubble))
            times.append(thubble)
        elif thubble - t_orig > 0.5*trange:
            if len(neigh_snaps) <= 2:
                # In the case that we need to extend a bit further
                neigh_snaps.append(ozy_name)
                tup = 0.5*(t_prev - thubble)
                weights.append(tup)
                times.append(thubble)
            break
    
    # And just go back to original place
    os.chdir(presentpath)

    neigh_snaps = np.asarray(neigh_snaps)
    weights = np.asarray(weights)

    if returnweight:
        return neigh_snaps,weights
    else:
        return neigh_snaps

def remove_out_zoom(obj, group):
    """Remove objects outside of zoom region.
    
    TODO: The details of the zoom region should be read from the simulation
            namelist.
    """
    # These details are for the NUT simulation.
    centre_zoom = obj.array(np.array([0.68,0.33,0.29]), 'code_length')
    radius_zoom = obj.quantity(0.34*0.5, 'code_length')
    
    add_group = True
    
    # Compute distance of group limits from zoom centre position.
    gdist_wrt_zoom = np.linalg.norm((group.position - centre_zoom))
    
    if gdist_wrt_zoom > radius_zoom:
        add_group = False
    
    return add_group

def info_printer(obj, grouptype, top):
    """General method to print data.

        TODO: Check and update group variables.
    """
    
    from ozy.group import grouptypes
    
    group_list = obj.__dict__[grouptypes[grouptype]]
    
    ngroups = len(group_list)
    
    if top > ngroups:
        top = ngroups
    
    time = 'z=%0.3f' % obj.simulation.redshift
    
    output  = '\n'
    output += '## Largest %d %s\n' % (top, grouptypes[grouptype])
    if hasattr(obj, 'data_file'): 
        output += '## from: %s\n' % obj.data_file
    output += '## %d @ %s' % (ngroups, time)
    output += '\n\n'
    
    cnt = 1
    
    if grouptype == 'halo':
        output += ' ID    Mdm       Mstar     Mgas      r         fgas\t|  CentralGalMstar\n'
        #         ' 0000  4.80e+09  4.80e+09  4.80e+09  7.64e-09  0.000\t|  7.64e-09'
        output += ' ---------------------------------------------------------------------------------\n'
        for o in group_list:
            cgsm = -1
            if (hasattr(o,'central_galaxy')) & (hasattr(o.central_galaxy,'masses')): 
                cgsm = o.central_galaxy.masses['stellar']
            output += ' %04d  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e\t|  %0.2e \n' % \
                      (o.ID, o.masses['dm'], o.masses['stellar'],
                       o.masses['gas'],o.radii['total_half_mass'],
                       cgsm)
            cnt += 1
            if cnt > top: 
                break
    elif grouptype == 'galaxy':
        output += ' ID    Mstar     Mgas      SFR       r         fgas   nrho      Central\t|  Mhalo     HID\n'
        output += ' ----------------------------------------------------------------------------------------\n'
        #         ' 0000  4.80e+09  4.80e+09  4.80e+09  7.64e-09  0.000  7.64e-09  False
        for o in group_list:
            phm, phid = -1, -1
            if o.halo is not None: phm, phid = o.halo.masses['total'], o.halo.GroupID
            output += ' %04d  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %s\t|  %0.2e  %d \n' % \
                      (o.GroupID, o.masses['stellar'], o.masses['gas'],
                       o.sfr, o.radii['total_half_mass'], 
                       o.local_number_density['1000'], o.central,
                       phm, phid)
            cnt += 1
            if cnt > top: break
    elif grouptype == 'cloud':
        output += ' ID    Mstar     Mgas      SFR       r         fgas   nrho      Central\t|  Mhalo     HID\n'
        output += ' ----------------------------------------------------------------------------------------\n'
        #         ' 0000  4.80e+09  4.80e+09  4.80e+09  7.64e-09  0.000  7.64e-09  False
        for o in group_list:
            halo = o.obj.galaxies[o.parent_galaxy_index].halo
            output += ' %04d  %0.2e  %0.2e  %0.2e  %0.2e   %0.2e  %s\t|  %0.2e  %d \n' % \
                      (o.GroupID, o.masses['stellar'], o.masses['gas'],
                       o.sfr, o.radii['total_half_mass'],
                       o.local_number_density['1000'], o.central,
                       halo.masses['dm'], halo.GroupID)
            cnt += 1
            if cnt > top: break

def RunningMedian(seq, M):
    """
     Purpose: Find the median for the points in a sliding window (odd number in size) 
              as it is moved from left to right by one point at a time.
      Inputs:
            seq -- list containing items for which a running median (in a sliding window) 
                   is to be calculated
              M -- number of items in window (window size) -- must be an integer > 1
      Otputs:
         medians -- list of medians with size N - M + 1
       Note:
         1. The median of a finite list of numbers is the "center" value when this list
            is sorted in ascending order. 
         2. If M is an even number the two elements in the window that
            are close to the center are averaged to give the median (this
            is not by definition)
    """   
    seq = iter(seq)
    s = []   
    m = M // 2

    # Set up list s (to be sorted) and load deque with first window of seq
    s = [item for item in islice(seq,M)]    
    d = deque(s)

    # Simple lambda function to handle even/odd window sizes    
    median = lambda : s[m] if bool(M&1) else (s[m-1]+s[m])*0.5

    # Sort it in increasing order and extract the median ("center" of the sorted window)
    s.sort()    
    medians = [median()]   

    # Now slide the window by one point to the right for each new position (each pass through 
    # the loop). Stop when the item in the right end of the deque contains the last item in seq
    for item in seq:
        old = d.popleft()          # pop oldest from left
        d.append(item)             # push newest in from right
        del s[bisect_left(s, old)] # locate insertion point and then remove old 
        insort(s, item)            # insert newest such that new sort is not required        
        medians.append(median())  
    return medians


def tidal_radius(central, satellite, method='BT87_simple'):
    """
    Computation of the tidal radius of a satellite with respect to a central galaxy.
    """
    from amr2 import amr_integrator,stats_utils
    from part2 import part_integrator

    if method == 'BT87_simple':
        """
        This calculation in Binney & Tremaine 1987 assumes that both the satellite and
        the host galaxy can be treated as singular isothermal spheres.
        
        The definition of tidal radius comes from approximating it as the Jacobi radius
        in Eq 8.107, and mutiplying by the factor 1.145 assuming that the host halo 
        can be approximated by a singular isothermal sphere.
        """
        # 1. Compute enclosed mass of the host halo within R0 (radial distance of satellite
        #    to the halo centre)

        output_path = central.obj.simulation.fullpath
        distance = central.position - satellite.position
        d = central.obj.quantity(np.linalg.norm(distance.to('code_length')),'code_length')

        # Initialise region
        selected_reg = init_region(central,'sphere',rmax=(d.to('kpc').d,'kpc'),rmin=(0.0,'kpc'))
        filt = init_filter('none','none',central)

        # Particles first
        glob_attrs = part_integrator.part_region_attrs()
        glob_attrs.nvars = 2
        glob_attrs.nwvars = 1
        part_integrator.allocate_part_regions_attrs(glob_attrs)
        glob_attrs.varnames.T.view('S128')[0] = b'dm/mass'.ljust(128)
        glob_attrs.varnames.T.view('S128')[1] = b'star/mass'.ljust(128)
        glob_attrs.wvarnames.T.view('S128')[0] = b'cumulative'.ljust(128)
        part_integrator.integrate_region(output_path,selected_reg,filt,glob_attrs)
        part_mass = central.obj.quantity(glob_attrs.data[0,0,0]+glob_attrs.data[1,0,0], 'code_mass')

        # Then gas
        glob_attrs = amr_integrator.amr_region_attrs()
        glob_attrs.nvars = 1
        glob_attrs.nfilter = 1
        amr_integrator.allocate_amr_regions_attrs(glob_attrs)
        glob_attrs.varnames.T.view('S128')[0] = b'mass'.ljust(128)
        glob_attrs.result[0].nbins = 5
        glob_attrs.result[0].nfilter = 1
        glob_attrs.result[0].nwvars = 1
        glob_attrs.result[0].varname = 'mass'
        mybins = get_code_bins(central.obj,'gas/mass',5)
        glob_attrs.result[0].scaletype = mybins[1]
        stats_utils.allocate_pdf(glob_attrs.result[0])
        glob_attrs.result[0].bins = mybins[0]
        glob_attrs.result[0].wvarnames.T.view('S128')[0] = b'mass'.ljust(128)
        glob_attrs.filters[0] = filt
        amr_integrator.integrate_region(output_path,selected_reg,False,glob_attrs)
        gas_mass = central.obj.quantity(glob_attrs.result[0].totweights[0,0], 'code_mass')

        tot_mass = gas_mass  + part_mass
        r = (satellite.virial_quantities['mass'] / (3.0*tot_mass))**(1.0/3.0) * d

    elif method == 'King62':
        #TODO: Add this calculation, which is significantly more complex
        pass
    else:
        print('This tidal radius method is not contemplated. Stoping!')
        exit
    
    return r

def structure_regions(group, position=None, radius=None,
                      add_substructure=True, add_neighbours=False,
                      add_all=False, add_intersections = False,
                      tidal_method='BT87_simple',rmax=(1e10,'kpc')):
    """
    This routine returns the regions of substructures so they can be used
    by the Ozymandias Fortran routines
    """
    from ozy.plot_settings import circle_dictionary
    from ozy.utils import tidal_radius
    
    mysubs = []

    # Get the central object information
    if group.type == 'halo':
        myhalo = group
    elif group.type == 'galaxy':
        myhalo = group.halo

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
                try:
                    tr  = s.radius[tidal_method]
                except:
                    tr = tidal_radius(myhalo,s,method=tidal_method)
                mysubs.append(init_region(mysub,'sphere',rmax=(tr.to('kpc'),'kpc'),
                            rmin=(0,'kpc')))
                
    # If asked for every structure in the halo finder, just add all
    if add_all:
        halos = group.obj.halos
        for h in halos:
            if h.npart >= 1000:
                r = h.virial_quantities['radius']
                mysubs.append(init_region(h,'sphere',rmax=(r.to('kpc'),'kpc'),
                            rmin=(0,'kpc')))
                
    # This looks for what virial spheres of other halos intersect with the
    # one provided. If position and radius are given, they're computed for
    # that instead of the group center and virial radius
    if add_intersections:
        if position == None and radius == None:
            position = group.position
            radius = rmax
        halos = group.obj.halos
        for h in halos:
            distance = position - h.position
            d = group.obj.quantity(np.linalg.norm(distance.to('kpc').d),'kpc')
            rsum = radius + h.virial_quantities['radius']
            if h.npart >= 1000 and rsum.to('kpc') >= d.to('kpc') and h.ID != myhalo.ID:
                r = h.virial_quantities['radius']
                mysubs.append(init_region(h,'sphere',rmax=(r.to('kpc'),'kpc'),
                            rmin=(0,'kpc')))
            
    # If asked for neighbours (so inside the virial radius) obtain them
    if add_neighbours:
        # TODO: This needs to be updated with the actual class procedure
        # myhalo.get_neighbours_in(rvir)
        pass
    return mysubs

def init_region(group, region_type, rmin=(0.0,'rvir'), rmax=(0.2,'rvir'), xmin=(0.0,'rvir'), xmax=(0.2,'rvir'),
                ymin=(0.0,'rvir'), ymax=(0.2,'rvir'),zmin=(0.0,'rvir'), zmax=(0.2,'rvir'),
                mycentre=([0.5,0.5,0.5],'rvir'), myaxis=np.array([1.,0.,0.]),
                return_enclosing_sphere=False):
    """Initialise region Fortran derived type with details of group."""
    from amr2 import vectors
    from amr2 import geometrical_regions as geo
    if not isinstance(rmin,tuple) or not isinstance(rmax,tuple):
        raise TypeError('The format for rmin and rmax should be %s, instead you gave for rmin %s and for rmax %s' %(type(tuple),type(rmin),type(rmax)))
        exit
    if not isinstance(zmin,tuple) or not isinstance(zmax,tuple):
        raise TypeError('The format for zmin and zmax should be %s, instead you gave for zmin %s and for zmax %s' %(type(tuple),type(zmin),type(zmax)))
        exit
    reg = geo.region()

    if region_type == 'sphere':
        reg.name = 'sphere'
        reg.criteria_name = 'r_sphere'
        centre = vectors.vector()
        centre.x, centre.y, centre.z = group.position[0], group.position[1], group.position[2]
        reg.centre = centre
        axis = vectors.vector()
        norm_L = group.angular_mom['total']/np.linalg.norm(group.angular_mom['total'])
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        reg.axis = axis
        bulk = vectors.vector()
        velocity = group.velocity.in_units('code_velocity')
        bulk.x, bulk.y, bulk.z = velocity[0].d, velocity[1].d, velocity[2].d
        reg.bulk_velocity = bulk
        if rmin[1] == 'rvir':
            reg.rmin = rmin[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.rmin = group.obj.quantity(rmin[0],str(rmin[1])).in_units('code_length')
        if rmax[1] == 'rvir':
            rmax = rmax[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            rmax = group.obj.quantity(rmax[0],str(rmax[1])).in_units('code_length')
        reg.rmax = rmax
        enclosing_sphere_p = group.position
        enclosing_sphere_r = rmax

    elif region_type == 'basic_sphere':
        reg.name = 'sphere'
        reg.criteria_name = 'r_sphere'
        centre = vectors.vector()
        centre.x, centre.y, centre.z = group.position[0], group.position[1], group.position[2]
        reg.centre = centre
        axis = vectors.vector()
        norm_L = group.angular_mom['total']/np.linalg.norm(group.angular_mom['total'])
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        reg.axis = axis
        bulk = vectors.vector()
        bulk.x, bulk.y, bulk.z = 0,0,0
        reg.bulk_velocity = bulk
        if rmin[1] == 'rvir':
            reg.rmin = rmin[0]*group.virial_quantities['radius'].d
        else:
            reg.rmin = group.obj.quantity(rmin[0],str(rmin[1])).in_units('code_length')
        if rmax[1] == 'rvir':
            rmax = rmax[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            rmax = group.obj.quantity(rmax[0],str(rmax[1])).in_units('code_length')
        reg.rmax = rmax
        enclosing_sphere_p = group.position
        enclosing_sphere_r = rmax
    elif region_type == 'basic_cube':
        reg.name = 'cube'
        centre = vectors.vector()
        centre.x, centre.y, centre.z = group.position[0], group.position[1], group.position[2]
        reg.centre = centre
        axis = vectors.vector()
        norm_L = group.angular_mom['total']/np.linalg.norm(group.angular_mom['total'])
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        reg.axis = axis
        bulk = vectors.vector()
        bulk.x, bulk.y, bulk.z = 0,0,0
        reg.bulk_velocity = bulk
        reg.xmin = group.obj.quantity(xmin[0],str(xmin[1])).in_units('code_length')
        reg.xmax = group.obj.quantity(xmax[0],str(xmax[1])).in_units('code_length')
        reg.ymin = group.obj.quantity(ymin[0],str(ymin[1])).in_units('code_length')
        reg.ymax = group.obj.quantity(ymax[0],str(ymax[1])).in_units('code_length')
        reg.zmin = group.obj.quantity(zmin[0],str(zmin[1])).in_units('code_length')
        reg.zmax = group.obj.quantity(zmax[0],str(zmax[1])).in_units('code_length')
        enclosing_sphere_p = group.position
        corner = np.array([reg.xmax-centre.x,reg.ymax-centre.y,reg.zmax-centre.z])
        enclosing_sphere_r = np.linalg.norm(corner)
        
    elif region_type == 'custom_cylinder':
        reg.name = 'cylinder'
        centre = vectors.vector()
        mycentre = group.obj.array(mycentre[0],str(mycentre[1])).in_units('code_length')
        centre.x, centre.y, centre.z = mycentre[0],mycentre[1],mycentre[2]
        reg.centre = centre
        axis = vectors.vector()
        norm_L = myaxis/np.linalg.norm(myaxis)
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        reg.axis = axis
        bulk = vectors.vector()
        velocity = group.velocity.in_units('code_velocity')
        bulk.x, bulk.y, bulk.z = velocity[0].d, velocity[1].d, velocity[2].d
        reg.bulk_velocity = bulk
        
        if rmin[1] == 'rvir':
            reg.rmin = rmin[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.rmin = group.obj.quantity(rmin[0],str(rmin[1])).in_units('code_length')
        if rmax[1] == 'rvir':
            reg.rmax = rmax[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.rmax = group.obj.quantity(rmax[0],str(rmax[1])).in_units('code_length')
        
        if zmin[1] == 'rvir':
            reg.zmin = zmin[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.zmin = group.obj.quantity(zmin[0],str(zmin[1])).in_units('code_length')
        if zmax[1] == 'rvir':
            reg.zmax = zmax[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.zmax = group.obj.quantity(zmax[0],str(zmax[1])).in_units('code_length')
        enclosing_sphere_p = mycentre + norm_L * max(rmax,0.5*(reg.zmax-reg.zmin))
        enclosing_sphere_r = np.sqrt(max(abs(reg.zmax),abs(reg.zmin))**2 + 2*reg.rmax**2)
    elif region_type == 'cylinder':
        reg.name = 'cylinder'
        centre = vectors.vector()
        centre.x, centre.y, centre.z = group.position[0], group.position[1], group.position[2]
        reg.centre = centre
        axis = vectors.vector()
        norm_L = group.angular_mom['total']/np.linalg.norm(group.angular_mom['total'])
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        reg.axis = axis
        bulk = vectors.vector()
        velocity = group.velocity.in_units('code_velocity')
        bulk.x, bulk.y, bulk.z = velocity[0].d, velocity[1].d, velocity[2].d
        reg.bulk_velocity = bulk

        if rmin[1] == 'rvir':
            reg.rmin = rmin[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.rmin = group.obj.quantity(rmin[0],str(rmin[1])).in_units('code_length')
        if rmax[1] == 'rvir':
            reg.rmax = rmax[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.rmax = group.obj.quantity(rmax[0],str(rmax[1])).in_units('code_length')
        
        if zmin[1] == 'rvir':
            reg.zmin = zmin[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.zmin = group.obj.quantity(zmin[0],str(zmin[1])).in_units('code_length')
        if zmax[1] == 'rvir':
            reg.zmax = zmax[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.zmax = group.obj.quantity(zmax[0],str(zmax[1])).in_units('code_length')
        enclosing_sphere_p = group.position + norm_L * max(rmax,0.5*(reg.zmax-reg.zmin))
        enclosing_sphere_r = np.sqrt(max(abs(reg.zmax),abs(reg.zmin))**2 + 2*reg.rmax**2)
    elif region_type == 'top_midplane_cylinder':
        reg.name = 'cylinder'
        axis = vectors.vector()
        norm_L = group.angular_mom['total']/np.linalg.norm(group.angular_mom['total'])
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        reg.axis = axis
        if rmin[1] == 'rvir':
            reg.rmin = rmin[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.rmin = group.obj.quantity(rmin[0],str(rmin[1])).in_units('code_length')
        if rmax[1] == 'rvir':
            reg.rmax = rmax[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.rmax = group.obj.quantity(rmax[0],str(rmax[1])).in_units('code_length')
        
        if zmin[1] == 'rvir':
            reg.zmin = zmin[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.zmin = group.obj.quantity(zmin[0],str(zmin[1])).in_units('code_length')
        if zmax[1] == 'rvir':
            reg.zmax = zmax[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.zmax = group.obj.quantity(zmax[0],str(zmax[1])).in_units('code_length')
        centre = vectors.vector()
        im_centre = group.position.in_units('code_length').d + norm_L.d * reg.zmax
        centre.x, centre.y, centre.z = im_centre[0], im_centre[1], im_centre[2]
        reg.centre = centre
        bulk = vectors.vector()
        velocity = group.velocity.in_units('code_velocity')
        bulk.x, bulk.y, bulk.z = velocity[0].d, velocity[1].d, velocity[2].d
        reg.bulk_velocity = bulk
        enclosing_sphere_p = im_centre
        enclosing_sphere_r = np.sqrt(max(abs(reg.zmax),abs(reg.zmin))**2 + 2*reg.rmax**2)
    elif region_type == 'bottom_midplane_cylinder':
        reg.name = 'cylinder'
        axis = vectors.vector()
        norm_L = -group.angular_mom['total']/np.linalg.norm(group.angular_mom['total'])
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        reg.axis = axis

        if rmin[1] == 'rvir':
            reg.rmin = rmin[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.rmin = group.obj.quantity(rmin[0],str(rmin[1])).in_units('code_length')
        if rmax[1] == 'rvir':
            reg.rmax = rmax[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.rmax = group.obj.quantity(rmax[0],str(rmax[1])).in_units('code_length')
        
        if zmin[1] == 'rvir':
            reg.zmin = zmin[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.zmin = group.obj.quantity(zmin[0],str(zmin[1])).in_units('code_length')
        if zmax[1] == 'rvir':
            reg.zmax = zmax[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d
        else:
            reg.zmax = group.obj.quantity(zmax[0],str(zmax[1])).in_units('code_length')
        centre = vectors.vector()
        im_centre = group.position.d + norm_L.d * reg.zmax
        centre.x, centre.y, centre.z = im_centre[0], im_centre[1], im_centre[2]
        reg.centre = centre
        bulk = vectors.vector()
        velocity = group.velocity.in_units('code_velocity')
        bulk.x, bulk.y, bulk.z = velocity[0].d, velocity[1].d, velocity[2].d
        reg.bulk_velocity = bulk
        enclosing_sphere_p = im_centre
        enclosing_sphere_r = np.sqrt(max(abs(reg.zmax),abs(reg.zmin))**2 + 2*reg.rmax**2)
    else:
        raise KeyError('Region type not supported. Please check!')
    if return_enclosing_sphere:
        return reg, enclosing_sphere_p, enclosing_sphere_r
    else:
        return reg

def init_filter(cond_strs, name, group):
    """Initialise filter Fortran derived type with the condition strings provided."""
    from amr2 import filtering

    if isinstance(cond_strs, str):
        cond_strs = [cond_strs]
    filt = filtering.filter()
    if cond_strs[0] == 'none' and name == 'none':
        filt.ncond = 0
        filt.name = 'none'
        return filt
    elif cond_strs[0] == 'none' and name != 'none':
        filt.ncond = 0
        filt.name = name
        return filt
    elif name != 'none':
        filt.ncond = len(cond_strs)
        filt.name = name
        filtering.allocate_filter(filt)
        for i in range(0, filt.ncond):
            # Variable name
            particle = False
            if cond_strs[i].split('/')[0].split('_')[0] == 'star' or cond_strs[i].split('/')[0].split('_')[0] == 'dm':
                correct_str = cond_strs[i].split('/')[0].split('_')[0] + '/' + '_'.join(cond_strs[i].split('/')[0].split('_')[1:])
                particle = True
            else:
                correct_str = cond_strs[i].split('/')[0]
            
            filt.cond_vars.T.view('S128')[i] = correct_str.ljust(128)
            # Expresion operator
            filt.cond_ops.T.view('S2')[i] = cond_strs[i].split('/')[1].ljust(2)
            # Value transformed to code units
            try:
                value = group.obj.quantity(float(cond_strs[i].split('/')[2]), cond_strs[i].split('/')[3])
                if particle:
                    filt.cond_vals[i] = value.in_units(get_code_units(correct_str.split('/')[1])).d
                else:
                    filt.cond_vals[i] = value.in_units(get_code_units(correct_str)).d
            except:
                # In the case of the condition value being a string
                # we use variables for the filters
                print('Using filter with variable instead of value!')
                filt.use_var[i] = True
                units1 = get_code_units(correct_str)
                units2 = get_code_units(cond_strs[i].split('/')[2])
                if units1 != units2:
                    raise ValueError("You cannot compare %s and %s"%(units1,units2))
                filt.cond_vars_comp.T.view('S128')[i] = cond_strs[i].split('/')[2].ljust(128)
                # And in place of units we should have the factor of that variable that we want
                filt.cond_vals[i] = cond_strs[i].split('/')[3]

        return filt
    else:
        raise ValueError("Condition strings are given, but not a name for the filter. Please set!")

def interp_nans(y):
    """
    This function fills up an array with NaNs by performing linear interpolation,
    useful when binning has become too small that some bins are empty.
    """
    nans, x = np.isnan(y), lambda z: z.nonzero()[0]
    y[nans]= np.interp(x(nans), x(~nans), y[~nans])

    nans, x = np.isinf(y), lambda z: z.nonzero()[0]
    y[nans]= np.interp(x(nans), x(~nans), y[~nans])

    return y

def get_SNevents_log(logfile,have_crs=False):
    """This routine allows a quick read of a RAMSES simulation output in which
        the MFB log has been activated.
    """
    import subprocess
    import os

    # Get creation time of file
    ts = os.path.getmtime(logfile)

    # Get line using grep
    result = subprocess.run('grep "MFB" '+str(logfile),shell=True,stdout=subprocess.PIPE)
    # Separate output in lines
    result = result.stdout.decode('utf-8').split('\n')

    if have_crs: 
        nvar = 8
        maxstr = 134
        print('Logfile with CRs!')
    else:
        nvar = 6
        maxstr = 74
    

    # Create data array
    sn_data = np.zeros((len(result)-1,nvar))

    if have_crs:

        for i in range(0, len(result)-1):
            line = result[i]
            if len(line) == maxstr:
                try:
                    sn_data[i,0] = float(line[6:16]) # z
                    sn_data[i,1] = float(line[20:28]) # N

                    sn_data[i,2] = float(line[42:50]) # nH
                    sn_data[i,3] = float(line[50:58]) # T
                    sn_data[i,4] = float(line[58:66]) # Z
                    sn_data[i,5] = float(line[66:74]) # eM
                    sn_data[i,6] = float(line[74:82]) # eCR

                    sn_data[i,7] = float(line[125:134]) # dx
                except:
                    print('Failure in line %i given by: '%i)
                    print(line)
    else:
        for i in range(0, len(result)-1):
            line = result[i]
            if len(line) == maxstr:
                try:
                    sn_data[i,0] = float(line[6:16]) # z
                    sn_data[i,1] = float(line[20:28]) # N

                    sn_data[i,2] = float(line[35:43]) # nH
                    sn_data[i,3] = float(line[43:51]) # T
                    sn_data[i,4] = float(line[51:59]) # Z

                    sn_data[i,5] = float(line[64:74]) # dx
                except:
                    print('Failure in line %i given by: '%i)
                    print(line)
    
    sn_data = sn_data[sn_data[:,0].argsort()[::-1]]

    print('Found %i SN events in %s'%(len(result)-1,logfile))

    return sn_data,ts

def sn_data_hdf5(logfiles,have_crs=False,outdir='Groups',filename='sn_catalogue.hdf5'):
    """
    This function builds the full catalogue of SN events found in a series of 
    RAMSES logfiles in which the MFB log is ON. Temporal cross-matching is done
    giving higher preference to more recent files.
    """
    import subprocess
    import os
    import h5py

    # Get variables depending on whether we have CRs or not
    if have_crs:
        variables = ['z','number','density','temperature','metallicity','magnetic_energy','cr_energy','dx']
    else:
        variables = ['z','number','density','temperature','metallicity','dx']

    # Check for the presence of a SN log catalogue
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    update_file = False
    if os.path.exists(os.path.join(outdir,filename)):
        print('SN catalogue file already present!')
        update_file = True
        snfile = h5py.File(os.path.join(outdir,filename),'r+')
        files_inside = snfile['files_used'].asstr()[:]
        # Only read the files that are not already included
        temp_logfiles = []
        for i in range(0, len(logfiles)):
            if not logfiles[i] in files_inside:
                temp_logfiles.append(logfiles[i])
        logfiles = temp_logfiles
        if len(logfiles) == 0:
            print('All asked files are already inside the found SN catalogue.')
            return snfile
    else:
        snfile = h5py.File(os.path.join(outdir,filename),'a')
        sim_name = os.getcwd().split('/')[-1]
        snfile.attrs.create('sim_name', sim_name)
        snfile.attrs.create('have_crs', have_crs)
    
    # Now load data from logfiles
    raw_sn_data = []
    log_timestamps = []
    files_read = []

    for i in range(0, len(logfiles)):
        if os.path.exists(logfiles[i]):
            sn_data, ts = get_SNevents_log(logfiles[i],have_crs=have_crs)
            if sn_data.shape[0] > 0:
                raw_sn_data.append(sn_data)
                log_timestamps.append(ts)
                files_read.append(str(logfiles[i]))
    if len(files_read) == 0:
            print('All asked files are already inside the found SN catalogue.')
            return snfile
    if not update_file:
        snfile.create_dataset('files_used',data=files_read,compression=1,maxshape=(None,))
    
    # Order files to figure out preference
    order_of_logs = np.asarray(log_timestamps).argsort()
    raw_sn_data = np.asarray(raw_sn_data)[order_of_logs]

    n_events = 0
    limits = np.zeros((len(raw_sn_data),2), dtype=int)

    for i in range(len(raw_sn_data),0,-1):
        if i==len(raw_sn_data):
            data = raw_sn_data[i-1]
            n_events += len(data)
            limits[i-1][0] = int(0)
            limits[i-1][1] = int(len(data))
        else:
            data = raw_sn_data[i-1]
            limits[i-1][0] = int(0)
            prev_limit = raw_sn_data[i][limits[i][0],0]
            z = data[:,0][::-1]
            limits[i-1][1] = int(z.searchsorted(prev_limit))
            n_events += limits[i-1][1]

    # Now, we can feed this into the full array
    full_sn_data = np.zeros((n_events,raw_sn_data[0].shape[1]))
    current_events = 0


    for i in range(0, len(raw_sn_data)):
        j,k = current_events,current_events+limits[i][1]
        data = raw_sn_data[i]
        full_sn_data[j:k,:] = data[limits[i][0]:limits[i][1],:]
        current_events += limits[i][1]

    # If there is data to add to an old file, make sure to not double count
    if update_file:
        orig_z = snfile['z'][:]
        low_limit = orig_z[::-1].searchsorted(full_sn_data[0,0])
        high_limit = orig_z[::-1].searchsorted(full_sn_data[-1,0])
        new_nevents = int(full_sn_data.shape[0] + high_limit + (len(orig_z)-low_limit))

        files_inside = snfile['files_used']
        new_files = np.concatenate((files_inside[:],files_read))
        files_inside.resize((len(new_files),))
        files_inside[:] = new_files

        for v in range(0, len(variables)):
            var = variables[v]
            orig_data = snfile[var]
            if (len(orig_z)-low_limit) == 0:
                new_data = np.concatenate((full_sn_data[:,v],orig_data[-high_limit:]))
            elif high_limit == 0:
                new_data = np.concatenate((orig_data[:(len(orig_z)-low_limit)],full_sn_data[:,v]))
            else:
                new_data = np.concatenate((orig_data[0:(len(orig_z)-low_limit)],full_sn_data[:,v]))
                new_data = np.concatenate((new_data, orig_data[:-high_limit]))
            orig_data.resize((new_nevents,))
            orig_data[:] = new_data

    else:
        for v in range(0, len(variables)):
            var = variables[v]
            data = full_sn_data[:,v]
            snfile.create_dataset(var,data=data,compression=1,maxshape=(None,))
    
    return snfile


def plot_cooling(cool_file):
    """
    This function allows for an easy inspection of the cooling curves
    saved in the RAMSES outputs and how they are used in post-processing.
    """
    from amr2 import cooling_module
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm,SymLogNorm
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from ozy.plot_settings import plotting_dictionary
    from ozy.dict_variables import common_variables,grid_variables,particle_variables,get_code_units
    
    # Read table
    cooling_module.read_cool(cool_file)
    mytable = cooling_module.cooling_table()
    cooling_module.retrieve_table(cool_file,mytable)

    Z = np.array([0.001,0.1,1.0])

    # Build figure
    figsize = plt.figaspect(float((6.0 * 1) / (5.0 * 3)))
    fig = plt.figure(figsize=figsize, facecolor='w',edgecolor='k')
    plot_grid = fig.add_gridspec(1, 3, wspace=0, hspace=0)
    ax = []
    for i in range(0,3):
        ax.append(fig.add_subplot(plot_grid[i]))
    ax = np.asarray(ax)

    # Solve cooling
    cooling_data = np.zeros((3,mytable.n1,mytable.n2))

    for i in range(0, mytable.n1):
        for j in range(0, mytable.n2):
            nH = 10**mytable.nh[i]
            T2 = 10**mytable.t2[j]
            for k in range(0,3):
                l,lp = cooling_module.solve_cooling(nH,T2,Z[k])
                cooling_data[k,i,j] = l*(nH**2)

    # Plot data
    for i in range(0,3):
        plotting_z = plotting_dictionary['net_cooling']
        ax[i].set_xlabel(r'$nH$ [cm$^{-3}$]',fontsize=18)
        if i==0:
            ax[i].set_ylabel(r'$T/\mu$ [K]',fontsize=18)
        else:
            ax[i].axes.yaxis.set_visible(False)
        
        ax[i].tick_params(labelsize=14,direction='in')
        ax[i].xaxis.set_ticks_position('both')
        ax[i].yaxis.set_ticks_position('both')
        ax[i].minorticks_on()
        ax[i].tick_params(which='major',axis="both",direction="in")

        ax[i].set_xscale('log')
        ax[i].set_yscale('log')
        x = 10**mytable.nh[:]
        y = 10**mytable.t2[:]
        z = cooling_data[i,:,:]
        print(abs(z).min(),z.max())
        plot = ax[i].pcolormesh(x,y,z.T,shading='auto',cmap=plotting_z['cmap'],
                                norm=SymLogNorm(linthresh=plotting_z['linthresh'],
                                linscale=plotting_z['linscale'],
                                vmin=plotting_z['vmin'],
                                vmax=plotting_z['vmax'],
                                base=10))
        ax[i].text(0.5, 0.9, r'$Z = %.3f Z_{\odot}$'%Z[i],
                            transform=ax[i].transAxes, fontsize=14,verticalalignment='top',
                            color='black')
        if i==0:
            cbaxes = inset_axes(ax[i], width="300%", height="5%", loc='upper left',
                                bbox_to_anchor=(0.0, 0., 1.0, 1.05),
                                bbox_transform=ax[i].transAxes,borderpad=0)
            cbar = fig.colorbar(plot, cax=cbaxes, orientation='horizontal')
            cbar.set_label(plotting_z['label'],fontsize=20)
            cbar.ax.tick_params(labelsize=10)
            cbaxes.xaxis.set_label_position('top')
            cbaxes.xaxis.set_ticks_position('top')
    fig.subplots_adjust(top=0.85,bottom=0.13,left=0.1,right=0.99)
    fig.savefig(cool_file.split('.out')[0]+'_table.png',format='png',dpi=300)

    return cooling_data


def gent_curve_rho(limit,T):
    from unyt import erg,g,K,cm
    s_hot = 23.2e+8 * erg / K / g
    s_cold = 4.4e+8 *  erg / K / g
    cv = 1.4e+8 * erg / K / g
    gamma = 5/3
    rho = 0.0
    if limit == 'hot':
        rho = 1.673532784796145e-24 * (T/(np.exp(s_hot/cv)*K)) ** (1/(gamma-1)) * g/cm**3
    elif limit == 'cold':
        rho = 1.673532784796145e-24 * (T/(np.exp(s_cold/cv)*K)) ** (1/(gamma-1)) * g/cm**3
    
    return rho

def gent_curve_T(limit,rho):
    from unyt import erg,g,K,cm
    s_hot = 23.2e+8 * erg / K / g
    s_cold = 4.4e+8 *  erg / K / g
    cv = 1.4e+8 * erg / K / g
    gamma = 5/3
    T = 0.0
    if limit == 'hot':
        T = (np.exp(s_hot/cv)*K) * (rho / (1.673532784796145e-24 * g/cm**3)) ** (gamma-1)
    elif limit == 'cold':
        T = (np.exp(s_cold/cv)*K) * (rho / (1.673532784796145e-24 * g/cm**3)) ** (gamma-1)
    
    return T


def stats_from_pdf(varname,x,PDF,xmin,xmax):
    """This function allows a quick computation of summary statistics
        used when normalised PDFs are returned from Ozymandias codes.
    """
    from scipy import interpolate, optimize

    # Check for misbehaving PDFs
    if any(np.isnan(PDF)):
        print('Empty PDF, ignoring!')
        return np.zeros(5)
    if all(PDF==0):
        print('Empty PDF, ignoring!')
        return np.zeros(5)
    if any(PDF<0):
        print('This PDF has negative values, so will be ignored!')
        # print(x,PDF)
        return np.zeros(5)
    if len(PDF[PDF!=0])==1:
        print('This PDF is composed of a single bin, so everything will be set to that value!')
        mean = x[PDF!=0][0]
        return np.array([mean,mean,0.0,mean,mean])

    CDF = np.cumsum(PDF)
    if CDF[-1] > 1.1 or CDF[-1]<0.9:
        print('Your PDF exceeds/lacks a total integral of 1 by more than 10%. Please check!')
        # print(x,PDF)
        PDF = PDF/np.sum(PDF)
        CDF = np.cumsum(PDF)
    
    # Mean using sum over x * PDF(x)
    mean = np.sum(x*PDF)

    # Standard deviation using sqrt of the sum over (x-mean)^2 * PDF
    std = np.sqrt(abs(np.sum((x-mean)**2*PDF)))
    # Get interpolation of discrete CDF
    f = interpolate.interp1d(x,CDF)
    
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots(1, 1, sharex=True, figsize=(6,4), dpi=100, facecolor='w', edgecolor='k')
    # ax.set_ylim([0,1])
    # ax.step(x,CDF)
    # ax.plot(x,f(x))
    # ax.plot([mean,mean],[0,1])
    # fig.savefig('/mnt/zfsusers/currodri/Codes/ozymandias/tests/test_cdf.png')
    # fig, ax = plt.subplots(1, 1, sharex=True, figsize=(6,4), dpi=100, facecolor='w', edgecolor='k')
    # ax.step(x,PDF)
    # ax.plot([mean,mean],[0,max(PDF)])
    # fig.savefig('/mnt/zfsusers/currodri/Codes/ozymandias/tests/test_pdf.png')
    interp_median = lambda x: f(x) - 0.5
    interp_q2 = lambda x: f(x) - 0.25
    interp_q4 = lambda x: f(x) - 0.75

    try:
        # Find roots of the interpolated CDF using the Newton-Rhapson method,
        # beginning from the mean value
        median = optimize.newton(interp_median,mean)
        q2 = optimize.newton(interp_q2,median)
        q4 = optimize.newton(interp_q4,median)
    except:
        # If that doesn't converge, just try the safer Brent's method
        print(interp_q2(min(x)),interp_q2(max(x)))
        if interp_q2(min(x)) >= 0.0:
            # In the case that everything fails...
            median = min(x)
            q2 = min(x)
            q4 = max(x[PDF>0.0])
        else:
            median = optimize.brentq(interp_median,min(x),max(x))
            q2 = optimize.brentq(interp_q2,min(x),max(x))
            q4 = optimize.brentq(interp_q4,min(x),max(x))
    if median > xmax or median < xmin:
        print('The median is out of bounds. Check!')
        print(varname)
        print(median,xmax,xmin)
        print(x[PDF!=0],PDF[PDF!=0])
        print(x,PDF)
        # exit(0)
    if q2 > xmax or q2 < xmin:
        print('Second quartile is out of bounds. Check!')
        print(varname)
        print(q2,xmax,xmin)
        print(x[PDF!=0],PDF[PDF!=0])
        # exit(0)
    if q4 > xmax or q4 < xmin:
        print('Fourth quartile is out of bounds. Check!')
        print(varname)
        print(q4,xmax,xmin)
        print(x[PDF!=0],PDF[PDF!=0])
        # exit(0)
    return np.array([mean,median,std,q2,q4])

def pdf_handler_to_stats(obj,pdf_obj,ifilt):
    from ozy.plot_settings import plotting_dictionary, \
                                symlog_variables
    # This returns:
    # mean,median,std,q2,q4,minvalue,max_value
    nwvar = pdf_obj.nwvars
    stats_array = np.zeros((nwvar,7))
    # Some fields have specific numerical flags at the end
    # which do not interfere with the units. If that is 
    # the case, get rid of that last 
    varname = str(pdf_obj.varname.decode("utf-8")).rstrip()
    scaletype = str(pdf_obj.scaletype.decode("utf-8")).rstrip() 
    plotting_def = plotting_dictionary[varname]
    try:
        numflag = int(varname.split('_')[-1])
        numflag = True
    except:
        numflag = False
    if numflag:
        sfrstr = varname.split('_')[0] +'_'+ varname.split('_')[1]
        code_units = get_code_units(sfrstr)
    else:
        code_units = get_code_units(varname)
    # print('ifilt',ifilt,pdf_obj.minv[:],pdf_obj.maxv[:],pdf_obj.minv[ifilt],pdf_obj.maxv[ifilt])
    for i in range(0, nwvar):
        PDF = pdf_obj.heights[ifilt,i,:]
        x = 0.5*(pdf_obj.bins[1:]+pdf_obj.bins[:-1])
        # print(varname,x,PDF,obj.array(np.array([pdf_obj.minv[ifilt],pdf_obj.maxv[ifilt]]),code_units))
        xmin,xmax = pdf_obj.minv[ifilt],pdf_obj.maxv[ifilt]
        if scaletype == 'log_even':
            xmin,xmax = np.log10(pdf_obj.minv[ifilt]),np.log10(pdf_obj.maxv[ifilt])
        stats_array[i,:5] = stats_from_pdf(varname,x,PDF,xmin,xmax)
        # print(varname,i,x,PDF)
        if not all(stats_array[i,:5]==0.0):
            # Just make sure no empty PDF
            if scaletype == 'log_even':
                # Propagation of errors from x to 10^x
                orig_sigma = stats_array[i,2]
                stats_array[i,:5] = 10**stats_array[i,:5]
                new_sigma = np.log(10)*orig_sigma*stats_array[i,0]
                stats_array[i,2] = new_sigma
            stats_array[i,5:] = np.array([pdf_obj.minv[ifilt],pdf_obj.maxv[ifilt]])
        elif any(np.isnan(PDF)) and (pdf_obj.minv[ifilt] != 0 or pdf_obj.maxv[ifilt] != 0):
            print('The limits of the binning may be wrong, because you have valid min and max!')
            # print(varname,obj.array(np.array([pdf_obj.minv[ifilt],pdf_obj.maxv[ifilt]]),code_units).to(plotting_def['units']),plotting_def['bin_min'],plotting_def['bin_max'])
            
    stats_array = obj.array(stats_array,code_units)
    return stats_array
    

def symlog_bins(min_val, max_val, n_bins, zero_eps=0.1, padding=0):
    """
    Splits a data range into log-like bins but with 0 and negative values taken into account.
    Can be used together with matplotlib 'symlog' axis sacale (i.e. ax.set_xscale('symlog'))
    Feel free to contribute: https://gist.github.com/artoby/0bcf790cfebed5805fbbb6a9853fe5d5
    """
    a = min_val / (1 + padding)
    b = max_val * (1 + padding)
        
    if a > b:
        a, b = b, a
        
    neg_range_log = None
    if a < -zero_eps:
        neg_range_log = [np.log10(-a), np.log10(zero_eps)]
    
    # Add a value to zero bin edges in case a lies within [-zero_eps; zero_eps) - so an additional bin will be added before positive range
    zero_bin_edges = []
    if -zero_eps <= a < zero_eps:
        zero_bin_edges = [a]
            
    pos_range_log = None
    if b > zero_eps:
        pos_range_log = [np.log10(max(a, zero_eps)), np.log10(b)]

    nonzero_n_bin_edges = n_bins + 1 - len(zero_bin_edges)
    
    neg_range_log_size = (neg_range_log[0] - neg_range_log[1]) if neg_range_log is not None else 0
    pos_range_log_size = (pos_range_log[1] - pos_range_log[0]) if pos_range_log is not None else 0
    
    range_log_size = neg_range_log_size + pos_range_log_size
    pos_n_bin_edges_raw = int(round(nonzero_n_bin_edges * (pos_range_log_size/range_log_size))) if range_log_size > 0 else 0
    # Ensure each range has at least 2 edges if it's not empty
    neg_n_bin_edges = max(2, nonzero_n_bin_edges - pos_n_bin_edges_raw) if neg_range_log_size > 0 else 0
    pos_n_bin_edges = max(2, nonzero_n_bin_edges - neg_n_bin_edges) if pos_range_log_size > 0 else 0
    
    neg_bin_edges = []
    if neg_n_bin_edges > 0:
        neg_bin_edges = list(-np.logspace(neg_range_log[0], neg_range_log[1], neg_n_bin_edges))
        
    pos_bin_edges = []
    if pos_n_bin_edges > 0:
        pos_bin_edges = list(np.logspace(pos_range_log[0], pos_range_log[1], pos_n_bin_edges))
    
    result = neg_bin_edges + zero_bin_edges + pos_bin_edges
    return np.asarray(result)

def get_code_bins(obj,varname,nbins=100,logscale=True):
    """This function provides bins for RAMSES variables in 
        code units, taking into account issues with variables
        with negative values."""
    from ozy.plot_settings import plotting_dictionary, \
                                symlog_variables
    from ozy.dict_variables import check_need_neighbours, common_variables, \
                                    grid_variables, \
                                    particle_variables, \
                                    get_code_units,basic_conv
    # Begin by checking the existence of the variable
    var_type = varname.split('/')[0]
    var_name = varname.split('/')[1]
    ok_var = False
    if var_type == 'gas':
        if var_name in common_variables or var_name in grid_variables:
                ok_var = True
        else:
            raise KeyError('This gas variable is not supported. Please check!: %s', varname)
    elif var_type == 'star':
        if var_name.split('_')[0] == 'sfr':
            if len(var_name.split('_')) == 3:
                sfr_name = var_name.split('_')[0] +'_'+var_name.split('_')[1]
            else:
                sfr_name = var_name.split('_')[0]
            if sfr_name in particle_variables:
                ok_var = True
            else:
                raise KeyError('This star variable is not supported. Please check!')
        else:
            if var_name in common_variables or var_name in particle_variables:
                ok_var = True
            else:
                raise KeyError('This star variable is not supported. Please check!')
    elif var_type == 'dm':
        if var_name in common_variables or var_name in particle_variables:
            ok_var = True
        else:
            raise KeyError('This DM variable is not supported. Please check!')
    
    if not ok_var:
        print('Your variable is not found!')
        exit

    # If everything is fine, we go and compute the bin edges
    if varname.split('/')[0] == 'star' or varname.split('/')[0] == 'dm':
        plotting_def = plotting_dictionary[varname.split('/')[0]+'_'+varname.split('/')[1]]
        stellar = True
    else:
        plotting_def = plotting_dictionary[varname.split('/')[1]]

    # The quantities should be in code units, so we transform them
    min_val = obj.quantity(plotting_def['bin_min'],plotting_def['units'])
    max_val = obj.quantity(plotting_def['bin_max'],plotting_def['units'])
    # Some fields have specific numerical flags at the end
    # which do not interfere with the units. If that is 
    # the case, get rid of that last bit
    try:
        numflag = int(varname.split('/')[1].split('_')[-1])
        numflag = True
    except:
        numflag = False
    if numflag:
        sfrstr = varname.split('/')[1].split('_')[0] +'_'+ varname.split('/')[1].split('_')[1]
        code_units = get_code_units(sfrstr)
    else:
        code_units = get_code_units(varname.split('/')[1])
    min_val = min_val.to(code_units).d
    max_val = max_val.to(code_units).d
    if logscale:
        if varname.split('/')[1] not in symlog_variables:
            bin_edges = np.linspace(np.log10(min_val),np.log10(max_val),nbins+1)
            scaletype = 'log_even'
        else:
            linscale = obj.quantity(plotting_def['linscale'],plotting_def['units'])
            bin_edges = symlog_bins(min_val,max_val,nbins,zero_eps=linscale.to(code_units).d)
            scaletype = 'symlog'
    else:
        bin_edges = np.linspace(min_val,max_val,nbins+1)
        scaletype = 'linear_even'
        
    return bin_edges,scaletype

    





    



