from matplotlib.pyplot import plot
import numpy as np
from collections import deque,Counter
from bisect import insort, bisect_left
from itertools import islice
import sys
import os
import subprocess
import re
import matplotlib.text as mtext
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
from unyt import unyt_array,unyt_quantity
from .variables_settings import geometrical_variables, raw_gas_variables, \
           raw_part_variables, derived_gas_variables, \
           derived_part_variables, star_variables, gravity_variables, \
           circle_dictionary, basic_conv, hydro_variables_ordering, \
           part_variables_ordering, part_variables_type

def ends_with_int_suffix(s):
    parts = s.rsplit('_', 1)
    if len(parts) == 2:
        return parts[1].isdigit()
    return False

def remove_last_suffix_if_numeric(s):
    if ends_with_int_suffix(s):
        return s.rsplit('_', 1)[0]
    return s

def get_code_units(varname,vartype):
    if vartype == 'gas':
        if varname in geometrical_variables:
            unit = geometrical_variables[varname]['code_units']
        elif varname in raw_gas_variables:
            unit = raw_gas_variables[varname]['code_units']
        elif varname in derived_gas_variables:
            unit = derived_gas_variables[varname]['code_units']
        elif varname in gravity_variables:
            unit = gravity_variables[varname]['code_units']
        else:
            raise KeyError('Gas variable not found, check: '+str(varname))
    elif vartype == 'part':
        if varname in geometrical_variables:
            unit = geometrical_variables[varname]['code_units']
        elif varname in raw_part_variables:
            unit = raw_part_variables[varname]['code_units']
        elif varname in derived_part_variables:
            unit = derived_part_variables[varname]['code_units']
        elif varname in star_variables:
            unit = star_variables[varname]['code_units']
        else:
            raise KeyError('Part variable not found, check: '+str(varname))
    return unit

def get_part_vartype(varname):
    if varname in part_variables_type:
        varttype = part_variables_type[varname]
    elif varname in geometrical_variables or \
        varname in raw_part_variables or \
        varname in derived_part_variables or \
        varname in star_variables:
        varttype = 1
    else:
        raise KeyError('Variable not found, check: '+str(varname))
    return varttype

    
def check_need_neighbours(varname,vartype):
    if vartype == 'gas':
        if varname in geometrical_variables:
            need = geometrical_variables[varname]['neighbour']
        elif varname in raw_gas_variables:
            need = raw_gas_variables[varname]['neighbour']
        elif varname in derived_gas_variables:
            need = derived_gas_variables[varname]['neighbour']
        elif varname in gravity_variables:
            need = gravity_variables[varname]['neighbour']
        else:
            raise KeyError('Gas variable not found, check: '+str(varname))
    elif vartype == 'part':
        if varname in geometrical_variables:
            need = geometrical_variables[varname]['neighbour']
        elif varname in raw_part_variables:
            need = raw_part_variables[varname]['neighbour']
        elif varname in derived_part_variables:
            need = derived_part_variables[varname]['neighbour']
        elif varname in star_variables:
            need = star_variables[varname]['neighbour']
        else:
            raise KeyError('Part variable not found, check: '+str(varname))
    return need

def check_need_gravity(varname,vartype):
    if vartype == 'gas':
        if varname in geometrical_variables:
            need = geometrical_variables[varname].get('gravity', False)
        elif varname in raw_gas_variables:
            need = raw_gas_variables[varname].get('gravity', False)
        elif varname in derived_gas_variables:
            need = derived_gas_variables[varname].get('gravity', False)
        elif varname in gravity_variables:
            need = gravity_variables[varname].get('gravity', False)
        else:
            raise KeyError('Gas variable not found, check: '+str(varname))
    elif vartype == 'part':
        if varname in geometrical_variables:
            need = geometrical_variables[varname].get('gravity', False)
        elif varname in raw_part_variables:
            need = raw_part_variables[varname].get('gravity', False)
        elif varname in derived_part_variables:
            need = derived_part_variables[varname].get('gravity', False)
        elif varname in star_variables:
            need = star_variables[varname].get('gravity', False)
        else:
            raise KeyError('Part variable not found, check: '+str(varname))
    return need

def check_need_rt(varname,vartype):
    if vartype == 'gas':
        if varname in geometrical_variables:
            need = geometrical_variables[varname].get('rt', False)
        elif varname in raw_gas_variables:
            need = raw_gas_variables[varname].get('rt', False)
        elif varname in derived_gas_variables:
            need = derived_gas_variables[varname].get('rt', False)
        elif varname in gravity_variables:
            need = gravity_variables[varname].get('rt', False)
        else:
            raise KeyError('Gas variable not found, check: '+str(varname))
    elif vartype == 'part':
        if varname in geometrical_variables:
            need = geometrical_variables[varname].get('rt', False)
        elif varname in raw_part_variables:
            need = raw_part_variables[varname].get('rt', False)
        elif varname in derived_part_variables:
            need = derived_part_variables[varname].get('rt', False)
        elif varname in star_variables:
            need = star_variables[varname].get('rt', False)
        else:
            raise KeyError('Part variable not found, check: '+str(varname))
    return need

def get_plotting_def(varname,vartype):
    if vartype == 'gas':
        if varname in geometrical_variables:
            plotting_def = geometrical_variables[varname]
        elif varname in raw_gas_variables:
            plotting_def = raw_gas_variables[varname]
        elif varname in derived_gas_variables:
            plotting_def = derived_gas_variables[varname]
        elif varname in gravity_variables:
            plotting_def = gravity_variables[varname]
        else:
            raise KeyError('Gas variable not found, check: '+str(varname))
    elif vartype == 'part':
        if varname in geometrical_variables:
            plotting_def = geometrical_variables[varname]
        elif varname in raw_part_variables:
            plotting_def = raw_part_variables[varname]
        elif varname in derived_part_variables:
            plotting_def = derived_part_variables[varname]
        elif varname in star_variables:
            plotting_def = star_variables[varname]
        else:
            raise KeyError('Part variable not found, check: '+str(varname))
    return plotting_def

def get_mu(X,Y):
    
    return 1./(2.*X + 3./4.*Y)

def get_electron_mu(X,Y):
    
    return 1./(X + 1./2.*Y)

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
    # If the luminance of the opposite color is greater, return the opposite color
    if opposite_luminance > luminance:
        return (opposite_red, opposite_green, opposite_blue, alpha)

    # Otherwise, return black or white depending on the luminance of the original color
    if luminance < 0.5:
        return (0, 0, 0, alpha) # Black
    else:
        return (1, 1, 1, alpha) # White


def invert_tick_colours(ax,var,vartype,type_scale,vmin=None,vmax=None,
                        linthresh=None,linscale=None,orientation='horizontal'):
    from matplotlib.colors import LogNorm,SymLogNorm
    from matplotlib import colormaps

    fig = plt.gcf()
    plotting_def = get_plotting_def(var,vartype)
    cmap = colormaps.get_cmap(plotting_def['cmap'])
    if vmin == None:
        vmin = plotting_def['vmin'+type_scale]
    if vmax == None:
        vmax = plotting_def['vmax'+type_scale]
    if linthresh == None and plotting_def['symlog']:
        linthresh = plotting_def['linthresh']
    if linscale == None and plotting_def['symlog']:
        linscale = plotting_def['linscale']
    if orientation == 'horizontal':
        ticks_pos = ax.get_xticks()
        ticks_labels = ax.get_xticklabels()
    else:
        ticks_pos = ax.get_yticks()
        ticks_labels = ax.get_yticklabels()
    if not plotting_def['symlog']:
        norm = LogNorm(vmin=vmin,
                         vmax=vmax,
                         clip=True)
        for tp,tl in zip(ticks_pos,ticks_labels):
            rgba = cmap(norm(10**tp))
            new_rgba = most_contrast_rgba(rgba)
            tl.set_color(new_rgba)
        fig.canvas.draw()
    else:
        norm = SymLogNorm(vmin=vmin,
                         vmax=vmax,
                         linthresh=linthresh,
                         linscale=linscale,
                         clip=True)
        for tp,tl in zip(ticks_pos,ticks_labels):
            rgba = cmap(norm(tp))
            new_rgba = most_contrast_rgba(rgba)
            tl.set_color(new_rgba)
        fig.canvas.draw()
        

def as_si(x, ndp):
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'${m:s}\times 10^{{{e:d}}}$'.format(m=m, e=int(e))
    
def read_infofile(infopath):
    info = {}
    info['aexp'] = 0.0
    info['redshift'] = 0.0
    info['omega_b'] = 0.0

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
            if (linehead == "boxlen")  : info['boxlen']  = float(newline[13:].strip())
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

def read_headerfile(filename):
    data = {}
    particle_fields = []
    particle_counts = {}
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    for i, line in enumerate(lines):
        # New style with "# Family Count" header
        if line.startswith("#") and "family" in line.lower() and "count" in line.lower():
            j = i + 1
            while j < len(lines) and not lines[j].lower().startswith("particle fields"):
                parts = lines[j].split()
                if len(parts) == 2:
                    family, count = parts
                    particle_counts[family.strip()] = int(count)
                j += 1
        # Old style with "Total number of ..." lines
        elif line.lower().startswith("total number of"):
            if "dark matter" in line.lower():
                key = "dark_matter_particles"
            elif "star" in line.lower():
                key = "star_particles"
            elif "sink" in line.lower():
                key = "sink_particles"
            elif "particles" in line.lower():
                key = "total_particles"
            else:
                continue
            try:
                count = int(lines[i + 1].strip())
                particle_counts[key] = count
            except (IndexError, ValueError):
                continue

        # Handle particle fields line
        if line.lower().startswith("particle fields"):
            if i + 1 < len(lines):
                particle_fields = lines[i + 1].split()

    # Assemble final output
    data['particle_counts'] = particle_counts
    data['particle_fields'] = particle_fields
    data['num_fields'] = len(particle_fields)
    return data


def read_rtinfo(fullpath):
    """Read RT descriptor file `info_rt_XXXXX.txt` for a snapshot directory.

    Parses the file and returns a Python dictionary with scalar fields and
    arrays similar in style to `read_infofile`.

    Returns `None` if the file is not present or the snapshot index cannot be
    determined.
    """
    try:
        index = int(fullpath[-5:])
    except Exception:
        return None

    rtfile = os.path.join(fullpath, 'info_rt_%05d.txt' % index)
    if not os.path.exists(rtfile):
        return None

    with open(rtfile, 'r') as fh:
        lines = [ln.rstrip('\n') for ln in fh]

    rt = {}
    # initialize defaults
    rt.update(dict(nrtvar=0, nions=0, ngroups=0, iions=0, rtdp=0,
                   x=0.0, y=0.0, scale_np=0.0, scale_pf=0.0,
                   rt_c_fraction=0.0, n_star=0.0, t2_star=0.0, g_star=0.0))

    # helper: map various possible headers to keys
    header_map = {
        'nrtvar': ('nrtvar', int),
        'nions': ('nions', int),
        'ngroups': ('ngroups', int),
        'iions': ('iions', int),
        'rtprecision': ('rtdp', int),
        'rtdp': ('rtdp', int),
        'x': ('x', float),
        'y': ('y', float),
        'unit_np': ('scale_np', float),
        'scale_np': ('scale_np', float),
        'unit_pf': ('scale_pf', float),
        'scale_pf': ('scale_pf', float),
        'rt_c_frac': ('rt_c_fraction', float),
        'rt_c_fraction': ('rt_c_fraction', float),
        'n_star': ('n_star', float),
        't2_star': ('t2_star', float),
        'g_star': ('g_star', float)
    }

    # First pass: read scalars from labelled lines by searching for key anywhere
    num_re = re.compile(r'[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?')
    int_re = re.compile(r'[+-]?\d+')
    for ln in lines[:200]:
        lnl = ln.lower()
        if len(lnl) < 1:
            continue
        for h, (key, caster) in header_map.items():
            if h in lnl:
                # try to find first numeric token after '=' or anywhere
                eqpos = ln.find('=')
                fragment = ln[eqpos+1:] if eqpos != -1 else ln
                token = None
                if caster is int:
                    m = int_re.search(fragment)
                    if m:
                        token = m.group(0)
                else:
                    m = num_re.search(fragment)
                    if m:
                        token = m.group(0)
                if token is not None:
                    try:
                        rt[key] = caster(token)
                    except Exception:
                        pass
                break

    nIons = int(rt.get('nions', 0))
    nGroups = int(rt.get('ngroups', 0))

    # helpers to collect numeric tokens across subsequent lines
    def tokens_from(start_idx):
        for j in range(start_idx, len(lines)):
            for tok in lines[j].split():
                yield tok

    # Find and parse groupl0, groupl1 and spec2group blocks
    float_re = num_re
    int_re = re.compile(r'[+-]?\d+')
    i = 0
    while i < len(lines):
        ln = lines[i]
        lnl = ln.lower()
        # groupl0 / groupl1: gather floats across possibly multiple lines
        if 'groupl0' in lnl or 'groupl1' in lnl:
            keyname = 'groupl0' if 'groupl0' in lnl else 'groupl1'
            want = nGroups if nGroups > 0 else 0
            vals = []
            # start from current line: take substring after '=' if present, else whole line
            j = i
            while len(vals) < want and j < len(lines):
                seg = lines[j]
                if '=' in seg:
                    seg = seg.split('=', 1)[1]
                for m in float_re.finditer(seg):
                    if len(vals) >= want:
                        break
                    try:
                        vals.append(float(m.group(0)))
                    except Exception:
                        vals.append(0.0)
                j += 1
            rt[keyname] = vals
            # continue scanning from where we left
            i = j
            continue

        # spec2group: integers mapping species to groups
        if 'spec2group' in lnl:
            want = nIons if nIons > 0 else 0
            vals = []
            j = i
            while len(vals) < want and j < len(lines):
                seg = lines[j]
                if '=' in seg:
                    seg = seg.split('=', 1)[1]
                for m in int_re.finditer(seg):
                    if len(vals) >= want:
                        break
                    try:
                        vals.append(int(m.group(0)))
                    except Exception:
                        vals.append(0)
                j += 1
            rt['spec2group'] = vals
            i = j
            continue

        i += 1

    # Parse per-group blocks
    rt['group_egy'] = []
    rt['group_csn'] = []
    rt['group_cse'] = []

    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.strip().lower().startswith('---group'):
            # begin parsing tokens after this line
            token_gen = tokens_from(i+1)
            # group energy: one value
            egy = 0.0
            try:
                t = next(token_gen)
                egy = float(t)
            except StopIteration:
                egy = 0.0
            except Exception:
                egy = 0.0

            # collect csn and cse: nIons each
            csn = []
            cse = []
            for _ in range(nIons):
                try:
                    t = next(token_gen)
                    csn.append(float(t))
                except StopIteration:
                    csn.append(0.0)
                except Exception:
                    csn.append(0.0)
            for _ in range(nIons):
                try:
                    t = next(token_gen)
                    cse.append(float(t))
                except StopIteration:
                    cse.append(0.0)
                except Exception:
                    cse.append(0.0)

            rt['group_egy'].append(egy)
            rt['group_csn'].append(csn)
            rt['group_cse'].append(cse)
            # advance i to continue after tokens we consumed: approximate by moving one line
        i += 1

    # Ensure keys exist
    if 'spec2group' not in rt:
        rt['spec2group'] = []
    if 'groupl0' not in rt:
        rt['groupl0'] = []
    if 'groupl1' not in rt:
        rt['groupl1'] = []

    return rt

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
    
def get_tdyn(galaxy,cgm=False):
        """
        Computes the dynamical time-scale tdyn as
        the time required for a test particle to complete
        one full orbit at 0.2 Rvir.

        tdyn = 2pi*sqrt(R^3/(GM))
        where we assume M = Mgas+Mstars*Mdm is measured 
        within 0.2 Rvir
        """
        from unyt import G
        
        if cgm:
            Mtot = galaxy.halo.virial_quantities['mass'] + galaxy.mass['baryon'] + galaxy.mass['halo_gas']
            r = 0.8*galaxy.obj.halos[galaxy.parent_halo_index].virial_quantities['radius']
            tdyn = 2*np.pi*np.sqrt(r**3/(G*Mtot))
        else:
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
    import random
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
        try:
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
            del sim
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
        except:
            print('Missing neighbour snapshot: ',ozy_name)
        
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
        del sim
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
    
    from .group import grouptypes
    
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
    from .integrators import integrate_hydro, integrate_part

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

        distance = central.position - satellite.position
        d = central.obj.quantity(np.linalg.norm(distance.to('code_length')),'code_length')

        # Particles first
        glob_attrs = integrate_part(central.obj,rmin=(0.0,'kpc'),
                                    rmax=(d.to('kpc').d,'kpc'),region_type='sphere',
                                    variables=['mass'],
                                    weights=['cumulative'],
                                    do_binning=[False])
        part_mass = central.obj.quantity(glob_attrs.result.total[0,0,0,0], 'code_mass')

        # Then gas
        glob_attrs = integrate_hydro(central.obj,rmin=(0.0,'kpc'),
                                    rmax=(d.to('kpc').d,'kpc'),region_type='sphere',
                                    variables=['mass'],
                                    weights=['cumulative'],
                                    do_binning=[False])
        gas_mass = central.obj.quantity(glob_attrs.result.total[0,0,0,0], 'code_mass')

        tot_mass = gas_mass  + part_mass
        r = (satellite.virial_quantities['mass'] / (3.0*tot_mass))**(1.0/3.0) * d

    elif method == 'King62':
        #TODO: Add this calculation, which is significantly more complex
        pass
    else:
        ValueError('Tidal radius method %s not implemented.' % method)
    return r

def structure_regions(group, position=None, radius=None,
                      add_substructure=True, add_neighbours=False,
                      add_all=False, add_intersections = False,
                      tidal_method='BT87_simple',rmax=(1e10,'kpc'),
                      verbose=False):
    """
    This routine returns the regions of substructures so they can be used
    by the Ozymandias Fortran routines
    """
    from .variables_settings import circle_dictionary
    from .utils import tidal_radius
    
    mysubs = []

    # Get the central object information
    if group.type == 'halo':
        myhalo = group
    elif group.type == 'galaxy':
        myhalo = group.halo

    if isinstance(rmax, tuple):
        rmax = group.obj.quantity(rmax[0],rmax[1])
    # If asked for substructure, obtain the substructure of the host halo
    subs_counter = 0
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
                subs_counter += 1
                
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
    inter_counter = 0
    if add_intersections:
        if isinstance(position,unyt_array):
            position = group.position
        if isinstance(rmax,unyt_quantity):
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
                inter_counter += 1
    if verbose: print('subs_counter,inter_counter:',subs_counter,inter_counter)
            
    # If asked for neighbours (so inside the virial radius) obtain them
    if add_neighbours:
        # TODO: This needs to be updated with the actual class procedure
        # myhalo.get_neighbours_in(rvir)
        pass
    return mysubs

def init_region(group, region_type, rmin=(0.0,'rvir'), rmax=(0.2,'rvir'), xmin=(0.0,'rvir'), xmax=(0.2,'rvir'),
                ymin=(0.0,'rvir'), ymax=(0.2,'rvir'),zmin=(0.0,'rvir'), zmax=(0.2,'rvir'),
                mycentre=([0.5,0.5,0.5],'rvir'), myaxis=np.array([0.,0.,0.]),
                return_enclosing_sphere=False):
    """Initialise region Fortran derived type with details of group."""
    from .amr.amr2_pkg import vectors
    from .amr.amr2_pkg import geometrical_regions as geo
    if not isinstance(rmin,tuple) or not isinstance(rmax,tuple):
        raise TypeError('The format for rmin and rmax should be %s, instead you gave for rmin %s and for rmax %s' %(type(tuple),type(rmin),type(rmax)))
        exit
    if not isinstance(zmin,tuple) or not isinstance(zmax,tuple):
        raise TypeError('The format for zmin and zmax should be %s, instead you gave for zmin %s and for zmax %s' %(type(tuple),type(zmin),type(zmax)))
        exit
    reg = geo.region()
    boxlen = 1 #group.obj.simulation.boxsize.to('code_length').d # In code units

    if region_type == 'sphere':
        reg.name = 'sphere'
        reg.criteria_name = 'r_sphere'
        centre = vectors.vector()
        centre.x, centre.y, centre.z = group.position[0].to('code_length')/boxlen, group.position[1].to('code_length')/boxlen, group.position[2].to('code_length')/boxlen
        reg.centre = centre
        axis = vectors.vector()
        norm_L = group.angular_mom['total']/np.linalg.norm(group.angular_mom['total'])
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        if not all(myaxis == 0.0):
            axis.x,axis.y,axis.z = myaxis[0], myaxis[1], myaxis[2]
        reg.axis = axis
        bulk = vectors.vector()
        velocity = group.velocity.in_units('code_velocity')
        bulk.x, bulk.y, bulk.z = velocity[0].d, velocity[1].d, velocity[2].d
        reg.bulk_velocity = bulk
        if rmin[1] == 'rvir':
            reg.rmin = rmin[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d/boxlen
        else:
            reg.rmin = group.obj.quantity(rmin[0],str(rmin[1])).in_units('code_length')/boxlen
        if rmax[1] == 'rvir':
            rmax = rmax[0]*group.obj.halos[group.parent_halo_index].virial_quantities['radius'].d/boxlen
        else:
            rmax = group.obj.quantity(rmax[0],str(rmax[1])).in_units('code_length')/boxlen
        reg.rmax = rmax
        enclosing_sphere_p = group.position.to('code_length')/boxlen
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
        centre.x, centre.y, centre.z = group.position[0]/boxlen, group.position[1]/boxlen, group.position[2]/boxlen
        reg.centre = centre
        axis = vectors.vector()
        norm_L = group.angular_mom['total']/np.linalg.norm(group.angular_mom['total'])
        axis.x,axis.y,axis.z = norm_L[0], norm_L[1], norm_L[2]
        reg.axis = axis
        bulk = vectors.vector()
        bulk.x, bulk.y, bulk.z = 0,0,0
        reg.bulk_velocity = bulk
        reg.xmin = group.obj.quantity(xmin[0],str(xmin[1])).in_units('code_length')/boxlen
        reg.xmax = group.obj.quantity(xmax[0],str(xmax[1])).in_units('code_length')/boxlen
        reg.ymin = group.obj.quantity(ymin[0],str(ymin[1])).in_units('code_length')/boxlen
        reg.ymax = group.obj.quantity(ymax[0],str(ymax[1])).in_units('code_length')/boxlen
        reg.zmin = group.obj.quantity(zmin[0],str(zmin[1])).in_units('code_length')/boxlen
        reg.zmax = group.obj.quantity(zmax[0],str(zmax[1])).in_units('code_length')/boxlen
        enclosing_sphere_p = group.position/boxlen
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
        enclosing_sphere_p = mycentre + group.obj.array(norm_L * max(reg.rmax,0.5*(reg.zmax-reg.zmin)),'code_length')
        enclosing_sphere_r = group.obj.quantity(np.sqrt(max(abs(reg.zmax),abs(reg.zmin))**2 + 2*reg.rmax**2),'code_length')
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

def init_filter_hydro(cond_strs, name, obj):
    """Initialise filter Fortran derived type with the condition strings provided."""
    from .amr.amr2_pkg import filtering
    if isinstance(cond_strs, str):
        cond_strs = [cond_strs]
    filt = filtering.filter_hydro()
    check_none = all(x == 'none' for x in cond_strs)
    ncond_real = len(cond_strs) - cond_strs.count('none')
    if check_none and name == 'none':
        filt.ncond = 0
        filt.name = 'none'
        filtering.allocate_filter_hydro(filt)
        return filt
    elif check_none and name != 'none':
        filt.ncond = 0
        filt.name = name
        filtering.allocate_filter_hydro(filt)
        return filt
    elif name != 'none':
        filt.ncond = ncond_real
        filt.name = name
        filtering.allocate_filter_hydro(filt)
        cond_strs = list(filter(lambda x: x != 'none', cond_strs))
        for i in range(0, filt.ncond):
            if cond_strs[i] != 'none':
                correct_str = cond_strs[i].split('/')[0]
                filt.cond_vars_name.T.view('S128')[i] = correct_str.ljust(128)
                # Expresion operator
                filt.cond_ops.T.view('S2')[i] = cond_strs[i].split('/')[1].ljust(2)
                # Value transformed to code units
                try:
                    value = obj.quantity(float(cond_strs[i].split('/')[2]), cond_strs[i].split('/')[3])
                    filt.cond_vals[i] = value.in_units(get_code_units(correct_str,'gas')).d
                except:
                    # In the case of the condition value being a string
                    # we use variables for the filters
                    filt.use_var[i] = True
                    units1 = get_code_units(correct_str,'gas')
                    units2 = get_code_units(cond_strs[i].split('/')[2],'gas')
                    if units1 != units2:
                        raise ValueError("You cannot compare %s and %s"%(units1,units2))
                    filt.cond_vars_comp_name = cond_strs[i].split('/')[2]
                    # And in place of units we should have the factor of that variable that we want
                    filt.cond_vals[i] = cond_strs[i].split('/')[3]

        return filt
    else:
        raise ValueError("Condition strings are given, but not a name for the filter. Please set!")

def init_filter_part(cond_strs, name, obj):
    """ Initialise filter_part Fortran derived type with the condition strings provided."""
    from part2 import filtering
    if isinstance(cond_strs, str):
        cond_strs = [cond_strs]
    filt = filtering.filter_part()
    check_none = all(x == 'none' for x in cond_strs)
    ncond_real = len(cond_strs) - cond_strs.count('none')
    if check_none and name == 'none':
        filt.ncond = 0
        filt.name = 'none'
        filtering.allocate_filter_part(filt)
        return filt
    elif check_none and name != 'none':
        filt.ncond = 0
        filt.name = name
        filtering.allocate_filter_part(filt)
        return filt
    elif name != 'none':
        filt.ncond = ncond_real
        filt.name = name
        filtering.allocate_filter_part(filt)
        cond_strs = list(filter(lambda x: x != 'none', cond_strs))
        for i in range(0, filt.ncond):
            if cond_strs[i] != 'none':
                correct_str = cond_strs[i].split('/')[0]
                if len(correct_str.split('_')) > 1:
                    nonum_str = correct_str.split('_')[0]
                else:
                    nonum_str = correct_str
                filt.cond_vars_name.T.view('S128')[i] = correct_str.ljust(128)
                # Expresion operator
                filt.cond_ops.T.view('S2')[i] = cond_strs[i].split('/')[1].ljust(2)
                # Value transformed to code units
                try:
                    value = obj.quantity(float(cond_strs[i].split('/')[2]), cond_strs[i].split('/')[3])
                    if get_part_vartype(nonum_str) == 1:
                        filt.cond_vals_d[i] = value.in_units(get_code_units(nonum_str,'part')).d
                    elif get_part_vartype(nonum_str) == 2:
                        filt.cond_vals_i[i] = int(value.in_units(get_code_units(nonum_str,'part')).d)
                    elif get_part_vartype(nonum_str) == 3:
                        filt.cond_vals_l[i] = int(value.in_units(get_code_units(nonum_str,'part')).d)
                except:
                    # In the case of the condition value being a string
                    # we use variables for the filters
                    filt.use_var[i] = True
                    units1 = get_code_units(nonum_str,'part')
                    units2 = get_code_units(nonum_str[i].split('/')[2],'part')
                    if units1 != units2:
                        raise ValueError("You cannot compare %s and %s"%(units1,units2))
                    filt.cond_vars_comp_name = cond_strs[i].split('/')[2]
                    # And in place of units we should have the factor of that variable that we want
                    if get_part_vartype(nonum_str) == 1:
                        filt.cond_vals_d[i] = float(cond_strs[i].split('/')[3])
                    elif get_part_vartype(nonum_str) == 2:
                        filt.cond_vals_i[i] = int(cond_strs[i].split('/')[3])
                    elif get_part_vartype(nonum_str) == 3:
                        filt.cond_vals_l[i] = int(cond_strs[i].split('/')[3])
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

def get_equilibrium_temp(cool_file,Z):
    
    from .amr.amr2_pkg import cooling_module
    
    # Read the table
    cooling_module.read_cool(cool_file)
    mytable = cooling_module.cooling_table()
    cooling_module.retrieve_table(cool_file,mytable)
    
    # Compute the net cooling
    cooling_data = np.zeros((mytable.n1,mytable.n2))
    for i in range(0, mytable.n1):
        for j in range(0, mytable.n2):
            nH = 10**mytable.nh[i]
            T2 = 10**mytable.t2[j]
            l,lp = cooling_module.solve_net_cooling(nH,T2,Z)
            cooling_data[i,j] = l*(nH**2)
    
    nH = 10.**mytable.nh[:]
    Tmu = 10.**mytable.t2[:]
    abs_lambda = abs(cooling_data.T)
    
    return nH,Tmu[np.argmin(abs_lambda,axis=0)]
    
    
def plot_cooling(cool_file):
    """
    This function allows for an easy inspection of the cooling curves
    saved in the RAMSES outputs and how they are used in post-processing.
    """
    from .amr.amr2_pkg import cooling_module
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm,SymLogNorm
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    
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
                l,lp = cooling_module.solve_net_cooling(nH,T2,Z[k])
                cooling_data[k,i,j] = l*(nH**2)

    # Plot data
    for i in range(0,3):
        plotting_z = get_plotting_def('net_cooling')
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
        # z = cooling_data[i,:,:]
        # print(abs(z).min(),z.max())
        # plot = ax[i].pcolormesh(x,y,z.T,shading='auto',cmap=plotting_z['cmap'],
        #                         norm=SymLogNorm(linthresh=plotting_z['linthresh'],
        #                         linscale=plotting_z['linscale'],
        #                         vmin=plotting_z['vmin'],
        #                         vmax=plotting_z['vmax'],
        #                         base=10))
        z = abs(cooling_data[i,:,:])
        print(abs(z).min(),z.max())
        plot = ax[i].pcolormesh(x,y,z.T,shading='auto',cmap=plotting_z['cmap'],
                                norm=LogNorm(vmin=1e-40,
                                vmax=1e-20))
        ax[i].text(0.5, 0.9, r'$Z = %.3f Z_{\odot}$'%Z[i],
                            transform=ax[i].transAxes, fontsize=14,verticalalignment='top',
                            color='black')
        absz = abs(z.T)
        ax[i].plot(x,y[np.argmin(absz,axis=0)],linestyle='--',color='k')
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
    from unyt import erg,g,K,cm,mp
    s_hot = 23.2e+8 * erg / K / g
    s_cold = 4.4e+8 *  erg / K / g
    cv = 1.4e+8 * erg / K / g
    gamma = 5/3
    rho = 0.0
    if limit == 'hot':
        rho = mp.to('g').d * (T/(np.exp(s_hot/cv)*K)) ** (1/(gamma-1)) * g/cm**3
    elif limit == 'cold':
        rho = mp.to('g').d * (T/(np.exp(s_cold/cv)*K)) ** (1/(gamma-1)) * g/cm**3
    
    return rho

def gent_curve_T(limit,rho):
    from unyt import erg,g,K,cm,mp
    s_hot = 23.2e+8 * erg / K / g
    s_cold = 4.4e+8 *  erg / K / g
    cv = 1.4e+8 * erg / K / g
    gamma = 5/3
    T = 0.0
    if limit == 'hot':
        T = (np.exp(s_hot/cv)*K) * (rho / (mp.to('g').d * g/cm**3)) ** (gamma-1)
    elif limit == 'cold':
        T = (np.exp(s_cold/cv)*K) * (rho / (mp.to('g').d * g/cm**3)) ** (gamma-1)
    
    return T


def stats_from_pdf(varname,x,PDF,xmin,xmax,verbose=False):
    """This function allows a quick computation of summary statistics
        used when normalised PDFs are returned from Ozymandias codes.
    """
    from scipy import interpolate, optimize

    # Check for misbehaving PDFs
    if any(np.isnan(PDF)):
        if verbose: print('Empty PDF, ignoring!')
        return np.zeros(5)
    if all(PDF==0):
        if verbose: print('Empty PDF, ignoring!')
        return np.zeros(5)
    if any(PDF<0):
        if verbose: print('This PDF has negative values, so will be ignored!')
        # print(x,PDF)
        return np.zeros(5)
    if len(PDF[PDF!=0])==1:
        if verbose: print('This PDF is composed of a single bin, so everything will be set to that value!')
        mean = x[PDF!=0][0]
        return np.array([mean,mean,0.0,mean,mean])

    CDF = np.cumsum(PDF)
    if CDF[-1] > 1.1 or CDF[-1]<0.9:
        if verbose: print('Your PDF exceeds/lacks a total integral of 1 by more than 10%. Please check!')
        # print(x,PDF)
        PDF = PDF/np.sum(PDF)
        CDF = np.cumsum(PDF)
    
    # Mean using sum over x * PDF(x)
    mean = np.sum(x*PDF)

    # Standard deviation using sqrt of the sum over (x-mean)^2 * PDF
    std = np.sqrt(abs(np.sum((x-mean)**2*PDF)))
    # Get interpolation of discrete CDF
    f = interpolate.interp1d(x,CDF)
    
    
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
        if verbose: print(interp_q2(min(x)),interp_q2(max(x)))
        if interp_q2(min(x)) >= 0.0:
            # In the case that everything fails...
            median = min(x)
            q2 = min(x)
            q4 = max(x[PDF>0.0])
        else:
            try:
                median = optimize.brentq(interp_median,min(x),max(x))
                q2 = optimize.brentq(interp_q2,min(x),max(x))
                q4 = optimize.brentq(interp_q4,min(x),max(x))
            except:
                print(varname,x,PDF,xmin,xmax)
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(1, 1, figsize=(6,4), dpi=100, facecolor='w', edgecolor='k')
                ax.set_ylim([0,1])
                ax.step(x,CDF)
                #ax.plot(x,f(x))
                #ax.plot([mean,mean],[0,1])
                fig.savefig('/mnt/zfsusers/currodri/Codes/ozymandias/tests/test_cdf.png')
                fig, ax = plt.subplots(1, 1, figsize=(6,4), dpi=100, facecolor='w', edgecolor='k')
                ax.step(x,PDF)
                ax.set_yscale('log')
                #ax.plot([mean,mean],[0,max(PDF)])
                fig.savefig('/mnt/zfsusers/currodri/Codes/ozymandias/tests/test_pdf.png')
                raise RuntimeError('All methods for the median failed. Please check! Figure saved in ozymandias/tests/test_pdf.png...')
    if median > xmax or median < xmin:
        if verbose: 
            print('The median is out of bounds. Check!')
            print(varname)
            print(median,xmax,xmin)
            print(x[PDF!=0],PDF[PDF!=0])
            print(x,PDF)
        # exit(0)
    if q2 > xmax or q2 < xmin:
        if verbose: 
            print('Second quartile is out of bounds. Check!')
            print(varname)
            print(q2,xmax,xmin)
            print(x[PDF!=0],PDF[PDF!=0])
        # exit(0)
    if q4 > xmax or q4 < xmin:
        if verbose: 
            print('Fourth quartile is out of bounds. Check!')
            print(varname)
            print(q4,xmax,xmin)
            print(x[PDF!=0],PDF[PDF!=0])
        # exit(0)
    return np.array([mean,median,std,q2,q4])

def pdf_handler_to_stats(obj,vartype,pdf_obj,ivar,ifilt,verbose=False):
    # This returns:
    # mean,median,std,q2,q4,minvalue,max_value
    nwvar = pdf_obj.nwvars
    stats_array = np.zeros((nwvar,7))
    # Some fields have specific numerical flags at the end
    # which do not interfere with the units. If that is 
    # the case, get rid of that last 
    varname = str(pdf_obj.varname.T.view('S128')[ivar][0].decode()).rstrip()
    scaletype = str(pdf_obj.scaletype.T.view('S128')[ivar][0].decode()).rstrip()
    clean_name = remove_last_suffix_if_numeric(varname)
    code_units = get_code_units(clean_name,vartype)
    for i in range(0, nwvar):
        wvarname = str(pdf_obj.wvarnames.T.view('S128')[i][0].decode()).rstrip()
        if wvarname != 'cumulative' and wvarname != 'counts':
            PDF = np.nan_to_num(pdf_obj.heights[ivar,ifilt,i,:],nan=0.0)
            x = 0.5*(pdf_obj.bins[1:,ivar]+pdf_obj.bins[:-1,ivar])
            if all(PDF==0.0):
                if pdf_obj.minv[ivar,ifilt] != 0 or pdf_obj.maxv[ivar,ifilt] != 0:
                    if verbose:
                        plt_def = get_plotting_def(varname,vartype)
                        minv = obj.quantity(pdf_obj.minv[ivar,ifilt],code_units)
                        maxv = obj.quantity(pdf_obj.maxv[ivar,ifilt],code_units)
                        minpdf = obj.quantity(pdf_obj.bins[0,ivar],code_units)
                        maxpdf = obj.quantity(pdf_obj.bins[-1,ivar],code_units)
                        print(f'Check the PDF limits ({minpdf.to(plt_def["units"])},{maxpdf.to(plt_def["units"])}), because you have valid min and max!: ',varname,minv.to(plt_def['units']),maxv.to(plt_def['units']))
                continue
            xmin,xmax = pdf_obj.minv[ivar,ifilt],pdf_obj.maxv[ivar,ifilt]
            if scaletype == 'log_even':
                if pdf_obj.minv[ivar,ifilt]<0.0 or pdf_obj.maxv[ivar,ifilt]<0.0:
                    print('Negative values in a log_even PDF: ',varname,pdf_obj.minv[ivar,ifilt],pdf_obj.maxv[ivar,ifilt])
                xmin,xmax = np.log10(pdf_obj.minv[ivar,ifilt]),np.log10(pdf_obj.maxv[ivar,ifilt])
            stats_array[i,:5] = stats_from_pdf(varname,x,PDF,xmin,xmax)
            if scaletype == 'log_even':
                # Propagation of errors from x to 10^x
                orig_sigma = stats_array[i,2]
                stats_array[i,:5] = 10**stats_array[i,:5]
                new_sigma = np.log(10)*orig_sigma*stats_array[i,0]
                stats_array[i,2] = new_sigma
            stats_array[i,5:] = np.array([pdf_obj.minv[ivar,ifilt],pdf_obj.maxv[ivar,ifilt]])
            if verbose and pdf_obj.nout[ivar,ifilt]/pdf_obj.nvalues[ivar,ifilt]>0.1: 
                print(f'Found points outside of PDF range for {varname}: {pdf_obj.nout[ivar,ifilt]}/{pdf_obj.nvalues[ivar,ifilt]}') 
                print(f'Range of values vs PDF limits: {xmin} - {xmax}, {x[0]} - {x[-1]}')
        else:
            stats_array[i,:] = pdf_obj.total[ivar,ifilt,i,0]
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
    return np.asarray(result),neg_n_bin_edges

def get_code_bins(obj,var_type,var_name,nbins=100,logscale=True,
                  minval=None,maxval=None,linthresh=None):
    """This function provides bins for RAMSES variables in 
        code units, taking into account issues with variables
        with negative values."""
    # Begin by checking the existence of the variable
    ok_var = False
    if var_type == 'gas':
        if var_name in geometrical_variables or var_name in raw_gas_variables \
            or var_name in derived_gas_variables or var_name in gravity_variables:
                ok_var = True
        else:
            raise KeyError('This gas variable is not supported. Please check!: %s', var_name)
    elif var_type == 'part':
        var_name_temp = remove_last_suffix_if_numeric(var_name)
        if var_name_temp in geometrical_variables or var_name_temp in raw_part_variables \
            or var_name_temp in derived_part_variables or var_name_temp in star_variables:
                ok_var = True
        else:
            raise KeyError('This part variable is not supported. Please check!: %s', var_name_temp)
    
    if not ok_var:
        print('Your variable is not found!')
        exit

    # If everything is fine, we go and compute the bin edges
    clean_name = remove_last_suffix_if_numeric(var_name)
    plotting_def = get_plotting_def(clean_name,var_type)

    # The quantities should be in code units, so we transform them
    if minval == None:
        min_val = obj.quantity(plotting_def['bin_min'],plotting_def['units'])
    else:
        min_val = minval
    if maxval == None:        
        max_val = obj.quantity(plotting_def['bin_max'],plotting_def['units'])
    else:
        max_val = maxval

    code_units = get_code_units(clean_name,var_type)
    min_val = min_val.to(code_units).d
    max_val = max_val.to(code_units).d
    zero_index = 0
    if logscale:
        if not plotting_def['symlog']:
            bin_edges = np.linspace(np.log10(min_val),np.log10(max_val),nbins+1)
            scaletype = 'log_even'
        else:
            if linthresh == None:
                linthresh = obj.quantity(plotting_def['linthresh'],plotting_def['units'])
            bin_edges,zero_index = symlog_bins(min_val,max_val,nbins,zero_eps=linthresh.to(code_units).d)
            scaletype = 'symlog'
            linthresh = linthresh.to(code_units).d
    else:
        bin_edges = np.linspace(min_val,max_val,nbins+1)
        scaletype = 'linear_even'
    if linthresh == None:
        linthresh = 0.0
    return bin_edges,scaletype,zero_index,linthresh