import numpy as np
from collections import deque,Counter
from bisect import insort, bisect_left
from itertools import islice
import sys
import os
import subprocess
from ozy.dict_variables import get_code_units

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

def closest_snap_z(simfolder,z):
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
        if t_orig - thubble <= 0.5*trange and t_orig > thubble:
            neigh_snaps.append(ozy_name)
            tup = 0.5*(t_orig - thubble)
            tdown = thubble - (t_orig - 0.5*trange)
            weights.append(tup+tdown)
            times.append(thubble)
        elif t_orig - thubble > 0.5*trange:
            if len(neigh_snaps) == 0:
                # In the case that we need to extend a bit further
                neigh_snaps.append(ozy_name)
                tup = 0.5*(t_orig - thubble)
                weights.append(tup)
                times.append(thubble)
                trange = 2*tup
            break
    
    # Add original snapshot
    neigh_snaps.append(ozy_orig)
    weights.append(0.0)
    times.append(t_orig)
    indexorig = len(weights)
    # And do the same for just above the original one
    for i in range(iorig+1,len(snapshots),1):
        ozy_name = 'ozy_%05d.hdf5'%(int(snapshots[i][-5:]))
        sim = ozy.load('Groups/'+ozy_name)
        cosmo = FlatLambdaCDM(H0=sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, 
                                        Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)
        thubble = cosmo.age(sim.simulation.redshift).value
        if thubble - t_orig <= 0.5*trange and t_orig < thubble:
            neigh_snaps.append(ozy_name)
            tdown = 0.5*(thubble - t_orig)
            tup = (t_orig + 0.5*trange) - thubble
            weights.append(tup+tdown)
            times.append(thubble)
        elif thubble - t_orig > 0.5*trange:
            if len(neigh_snaps) <= 2:
                # In the case that we need to extend a bit further
                neigh_snaps.append(ozy_name)
                tup = 0.5*(thubble - t_orig)
                weights.append(tup)
                times.append(thubble)
            break
    # And compute the weight of the original/middle snapshot
    tup = 0.5*(times[indexorig] - times[indexorig-1])
    tdown = 0.5*(times[indexorig-1] - times[indexorig-2])
    weights[indexorig-1] = tup+tdown
    
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

    if method == 'BT87_simple':
        r = (satellite.virial_quantities['mass']/(2.0*central.virial_quantities['mass']))**(1.0/3.0)
        r = r * satellite.virial_quantities['radius']

    elif method == 'BT87_centrifugal':
        mass_ratio = satellite.virial_quantities['mass']/central.virial_quantities['mass']
        r = (mass_ratio/(3.0 + mass_ratio))**(1.0/3.0)
        r = r * satellite.virial_quantities['radius']
    elif method == 'King62':
        #TODO: Add this calculation, which is significantly more complex
        pass
    else:
        print('This tidal radius method is not contemplated. Stoping!')
        exit
    
    return r

def init_region(group, region_type, rmin=(0.0,'rvir'), rmax=(0.2,'rvir'), xmin=(0.0,'rvir'), xmax=(0.2,'rvir'),
                ymin=(0.0,'rvir'), ymax=(0.2,'rvir'),zmin=(0.0,'rvir'), zmax=(0.2,'rvir')):
    """Initialise region Fortran derived type with details of group."""
    sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/amr')
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

    elif region_type == 'basic_sphere':
        reg.name = 'sphere'
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
            reg.rmax = rmax[0]*group.virial_quantities['radius'].d
        else:
            reg.rmax = group.obj.quantity(rmax[0],str(rmax[1])).in_units('code_length')
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
        im_centre = group.position.in_units('code_length').d + 0.99 * norm_L.d * reg.zmax
        centre.x, centre.y, centre.z = im_centre[0], im_centre[1], im_centre[2]
        reg.centre = centre
        bulk = vectors.vector()
        velocity = group.velocity.in_units('code_velocity')
        bulk.x, bulk.y, bulk.z = velocity[0].d, velocity[1].d, velocity[2].d
        reg.bulk_velocity = bulk
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
        im_centre = group.position.d + 0.99 * norm_L.d * reg.zmax
        centre.x, centre.y, centre.z = im_centre[0], im_centre[1], im_centre[2]
        reg.centre = centre
        bulk = vectors.vector()
        velocity = group.velocity.in_units('code_velocity')
        bulk.x, bulk.y, bulk.z = velocity[0].d, velocity[1].d, velocity[2].d
        reg.bulk_velocity = bulk
    else:
        raise KeyError('Region type not supported. Please check!')
    return reg

def init_filter(cond_strs, name, group):
    """Initialise filter Fortran derived type with the condition strings provided."""
    sys.path.append('/mnt/zfsusers/currodri/Codes/ozymandias/ozy/amr')
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
            filt.cond_vars.T.view('S128')[i] = cond_strs[i].split('/')[0].ljust(128)
            # Expresion operator
            filt.cond_ops.T.view('S2')[i] = cond_strs[i].split('/')[1].ljust(2)
            # Value transformed to code units
            value = group.obj.quantity(float(cond_strs[i].split('/')[2]), cond_strs[i].split('/')[3])
            filt.cond_vals[i] = value.in_units(get_code_units(cond_strs[i].split('/')[0])).d
        return filt
    else:
        raise ValueError("Condition strings are given, but a name for the filter. Please set!")

