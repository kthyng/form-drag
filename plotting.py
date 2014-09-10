'''
Plots from Admiralty Inlet 65 meter simulation.
'''

import scipy.io
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import visclaw.colormaps as colormaps
from skimage import color
from matplotlib import cm
import matplotlib as mpl

mpl.rcParams.update({'font.size': 14})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'


### Plot Bathymetry ###

# download bathymetry, which can be found at: http://figshare.com/preview/_preview/1165560 (27.3MB)

# Read in bathymetry
mat = scipy.io.loadmat('cascadia_gridded.mat')

# x and y limits for this plot
lonlims = [-122.8, -122.55]
latlims = [47.9665, 48.227]

# Functionality copied from https://github.com/clawpack/geoclaw/blob/master/src/python/geoclaw/topotools.py#L873
land_cmap = plt.get_cmap('Greens_r')
sea_cmap = plt.get_cmap('Blues_r')
cmap = colormaps.add_colormaps((land_cmap, sea_cmap), 
                                data_limits=[-200,175],
                                data_break=0.0)

# # test lightness profile of my colormap
# x = np.linspace(0.0, 1.0, 100) # indices to step through colormap
# # Get rgb values for colormap
# rgb = cmap(x)[np.newaxis,:,:3]

# # Get colormap in CIE LAB. We want the L here.
# lab = color.rgb2lab(rgb)

# fig = plt.figure(figsize=(11.5,4))
# ax = fig.add_subplot(111)
# ax.scatter(x, lab[0,::-1,0], c=x, cmap=cmap, s=300, linewidths=0.)
# ax.axis([-0.1,4.7,0,100])

# levels to plot
levs = np.concatenate((np.arange(-200, 0, 20), np.arange(0,200,25)))

# Make plot
fig = plt.figure(figsize=(9,8))
ax = fig.add_subplot(111)
mappable = ax.contourf(mat['lon_topo'], mat['lat_topo'], mat['z_topo'], cmap=cmap, levels=levs)
ax.set_xlim(lonlims)
ax.set_ylim(latlims)
ax.set_xlabel('Longitude [degrees]')
ax.set_ylabel('Latitude [degrees]')
# Turn off annoying offset, from https://github.com/clawpack/geoclaw/blob/master/src/python/geoclaw/topotools.py#L844
ax.ticklabel_format(format="plain", useOffset=False)
plt.xticks(rotation=20)
cb = fig.colorbar(mappable)
cb.set_label('Height/depth [m]')
plt.tight_layout()

# Save figure
fig.savefig('figures/domain.png')
fig.show()