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
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

mpl.rcParams.update({'font.size': 16})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'


### Plot Bathymetry of Puget Sound, Admiralty Inlet, and Admiralty Head ###

# download bathymetry, which can be found at: http://figshare.com/preview/_preview/1165560 (27.3MB)

# Read in bathymetry
mat = scipy.io.loadmat('cascadia_gridded.mat')

# x and y limits for these plots
lonlimsPS = [-123.21, -122.15];
latlimsPS = [47.02, 48.82];
lonlimsAI = [-122.8, -122.54]
latlimsAI = [47.9665, 48.227]
lonlimsAH = [-122.71, -122.65]
latlimsAH = [48.12, 48.18]

# Functionality copied from https://github.com/clawpack/geoclaw/blob/master/src/python/geoclaw/topotools.py#L873
land_cmap = plt.get_cmap('Greens_r')
sea_cmap = plt.get_cmap('Blues_r')
cmapPS = colormaps.add_colormaps((land_cmap, sea_cmap), 
                                data_limits=[-325,2500],
                                data_break=0.0)
cmapAI = colormaps.add_colormaps((land_cmap, sea_cmap), 
                                data_limits=[-200,175],
                                data_break=0.0)
cmapAH = colormaps.add_colormaps((land_cmap, sea_cmap), 
                                data_limits=[-110,50],
                                data_break=0.0)

# levels to plot
levsPS = np.concatenate((np.arange(-325, 0, 25), np.arange(0,3000,500)))
levsAI = np.concatenate((np.arange(-200, 0, 20), np.arange(0,350,175))) #200,25)))
levsAH = np.concatenate((np.arange(-120, 0, 20), np.arange(0,100,50)))


# Make Puget Sound plot
fig = plt.figure(figsize=(16,16))
axPS = fig.add_subplot(111)
mappablePS = axPS.contourf(mat['lon_topo'], mat['lat_topo'], mat['z_topo'], cmap=cmapPS, levels=levsPS)
# outline coast in case plot is printed
axPS.contour(mat['lon_topo'], mat['lat_topo'], mat['z_topo'], [0], lw=3, colors='0.15')
axPS.set_xlim(lonlimsPS)
axPS.set_ylim(latlimsPS)
axPS.set_xlabel('Longitude [degrees]')
axPS.set_ylabel('Latitude [degrees]')
# Turn off annoying offset, from https://github.com/clawpack/geoclaw/blob/master/src/python/geoclaw/topotools.py#L844
axPS.ticklabel_format(format="plain", useOffset=False)
plt.xticks(rotation=20)
cbPS = fig.colorbar(mappablePS)
cbPS.set_label('Height/depth [m]')
plt.tight_layout()
# Label
axPS.text(0.7, 0.025, 'Puget Sound', transform=axPS.transAxes, color='0.15')


# Inset magnified plot of Admiralty Inlet
axAI = zoomed_inset_axes(axPS, 2, loc=1)
mappableAI = axAI.contourf(mat['lon_topo'], mat['lat_topo'], mat['z_topo'], cmap=cmapAI, levels=levsAI)
# outline coast in case plot is printed
axAI.contour(mat['lon_topo'], mat['lat_topo'], mat['z_topo'], [0], lw=3, colors='0.15')
axAI.set_xlim(lonlimsAI)
axAI.set_ylim(latlimsAI)
# turn off ticks
plt.xticks(visible=False)
plt.yticks(visible=False)
plt.setp(axAI,xticks=[],yticks=[])
# Inlaid colorbar
caxAI = fig.add_axes([0.735, 0.77, 0.0125, 0.2])
cbAI = plt.colorbar(mappableAI, cax=caxAI, orientation='vertical')
# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(axPS, axAI, loc1=2, loc2=4, fc="none", ec="0.3", lw=1.5)
# Label
axAI.text(0.044, 0.04, 'Admiralty Inlet', transform=axAI.transAxes, color='0.15')


# Inset magnified plot of Admiralty Head
axAH = zoomed_inset_axes(axPS, 8, loc=3)
mappableAH = axAH.contourf(mat['lon_topo'], mat['lat_topo'], mat['z_topo'], cmap=cmapAH, levels=levsAH)
# outline coast in case plot is printed
axAH.contour(mat['lon_topo'], mat['lat_topo'], mat['z_topo'], [0], lw=3, colors='0.15')
axAH.set_xlim(lonlimsAH)
axAH.set_ylim(latlimsAH)
# turn off ticks
plt.xticks(visible=False)
plt.yticks(visible=False)
plt.setp(axAH,xticks=[],yticks=[])
# Inlaid colorbar
caxAH = fig.add_axes([0.35, 0.1, 0.0125, 0.2])
cbAH = plt.colorbar(mappableAH, cax=caxAH, orientation='vertical')
# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(axPS, axAH, loc1=2, loc2=4, fc="none", ec="0.3", lw=1.5)
# Label
axAH.text(0.45, 0.92, 'Admiralty Head', transform=axAH.transAxes, color='0.15')

plt.draw()
plt.show()

# Save figure
fig.savefig('figures/domains.png')
fig.show()


# ### Plot Bathymetry of just Admiralty Inlet ###

# # download bathymetry, which can be found at: http://figshare.com/preview/_preview/1165560 (27.3MB)

# # Read in bathymetry
# mat = scipy.io.loadmat('cascadia_gridded.mat')

# # x and y limits for this plot
# lonlims = [-122.8, -122.55]
# latlims = [47.9665, 48.227]

# # Functionality copied from https://github.com/clawpack/geoclaw/blob/master/src/python/geoclaw/topotools.py#L873
# land_cmap = plt.get_cmap('Greens_r')
# sea_cmap = plt.get_cmap('Blues_r')
# cmap = colormaps.add_colormaps((land_cmap, sea_cmap), 
#                                 data_limits=[-200,175],
#                                 data_break=0.0)

# # # test lightness profile of my colormap
# # x = np.linspace(0.0, 1.0, 100) # indices to step through colormap
# # # Get rgb values for colormap
# # rgb = cmap(x)[np.newaxis,:,:3]

# # # Get colormap in CIE LAB. We want the L here.
# # lab = color.rgb2lab(rgb)

# # fig = plt.figure(figsize=(11.5,4))
# # ax = fig.add_subplot(111)
# # ax.scatter(x, lab[0,::-1,0], c=x, cmap=cmap, s=300, linewidths=0.)
# # ax.axis([-0.1,4.7,0,100])

# # levels to plot
# levs = np.concatenate((np.arange(-200, 0, 20), np.arange(0,200,25)))

# # Make plot
# fig = plt.figure(figsize=(9,8))
# ax = fig.add_subplot(111)
# mappable = ax.contourf(mat['lon_topo'], mat['lat_topo'], mat['z_topo'], cmap=cmap, levels=levs)
# ax.set_xlim(lonlims)
# ax.set_ylim(latlims)
# ax.set_xlabel('Longitude [degrees]')
# ax.set_ylabel('Latitude [degrees]')
# # Turn off annoying offset, from https://github.com/clawpack/geoclaw/blob/master/src/python/geoclaw/topotools.py#L844
# ax.ticklabel_format(format="plain", useOffset=False)
# plt.xticks(rotation=20)
# cb = fig.colorbar(mappable)
# cb.set_label('Height/depth [m]')
# plt.tight_layout()

# # Save figure
# fig.savefig('figures/domain.png')
# f