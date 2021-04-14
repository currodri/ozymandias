import numpy as np
import sys
sys.path.append('/mnt/extraspace/currodri/Codes/ozymandias/ozy/amr')
import amr2 as amr2prof
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
sns.set(style="white")
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
hfont = {'fontname':'Helvetica'}
matplotlib.rc('text', usetex = True)
matplotlib.rc('font', **{'family' : "serif"})
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
matplotlib.rcParams.update(params)

repository = '/mnt/extraspace/currodri/NUT/cosmoNUThd/output_00044'
reg = amr2prof.geometrical_regions.region()
reg.name = 'cylinder'
centre = amr2prof.vectors.vector()
centre.x,centre.y,centre.z = 0.70336995, 0.33568012, 0.341516
reg.centre = centre
axis = amr2prof.vectors.vector()
axis.x,axis.y,axis.z = -0.35676077,  0.40252785,  0.84302614
reg.axis = axis
bulk = amr2prof.vectors.vector()
bulk.x,bulk.y,bulk.z = 28.26907349, -20.34032822, -25.46469116
reg.bulk_velocity = bulk
reg.rmin = 0.0
reg.rmax = 0.0028299868106842
reg.zmin = -0.0028299868106842
reg.zmax = 0.0028299868106842

filt = amr2prof.filtering.filter()
filt.ncond = 1
filt.name = 'outflow'
amr2prof.filtering.allocate_filter(filt)
filt.cond_vars.T.view('S128')[0] = b'v_cyl_r'.ljust(128)
filt.cond_ops.T.view('S2')[0] = b'>'.ljust(2)
filt.cond_vals[0] = 0.00250019

print(filt.cond_vars.T.view('S128')[0])

prof_data = amr2prof.amr_profiles.profile_handler()
prof_data.profdim=2
prof_data.xvarname='r_cyl'
prof_data.nyvar=3
prof_data.nwvar=2
prof_data.nbins=100
amr2prof.amr_profiles.allocate_profile_handler(prof_data)
prof_data.yvarnames.T.view('S128')[0] = b'density'.ljust(128)
prof_data.yvarnames.T.view('S128')[1] = b'v_cyl_phi'.ljust(128)
prof_data.yvarnames.T.view('S128')[2] = b'thermal_energy'.ljust(128)
ynames_plot = [r'$\rho$ [code units]',r'$v_{\phi}^{\rm cyl}$ [code units]', r'$E_{\rm th}$ [code units]']
prof_data.wvarnames.T.view('S128')[0] = b'mass'.ljust(128)
prof_data.wvarnames.T.view('S128')[1] = b'volume'.ljust(128)
amr2prof.amr_profiles.onedprofile(repository,reg,filt,prof_data,0)

# print(prof_data.ydata)
fig, axes = plt.subplots(prof_data.ydata.shape[1], 1, sharex=True, figsize=(5,7), dpi=100, facecolor='w', edgecolor='k')
for i in range(0, prof_data.ydata.shape[1]):
    ax = axes[i]
    ax.set_xlabel(r'$r$ [codelength]', fontsize=16)
    ax.set_ylabel(ynames_plot[i], fontsize=12)
    ax.set_xscale('log')
    ax.tick_params(labelsize=12)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.minorticks_on()
    ax.tick_params(which='both',axis="both",direction="in")
    ax.plot(prof_data.xdata, prof_data.ydata[:,i,0,0], marker='o', markersize=2)
fig.subplots_adjust(top=0.91, bottom=0.08,right=0.97,hspace=0.0,left=0.2)
fig.savefig('profile_NUT.png', format='png', dpi=200)

