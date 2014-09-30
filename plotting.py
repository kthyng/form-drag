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
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import MaxNLocator

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
lonlimsPS = np.array([-124., -122.15]) #-123.21, -122.15])
latlimsPS = np.array([47.02, 48.82])
lonlimsAI = np.array([-122.85, -122.535])
latlimsAI = np.array([47.9665, 48.228])
lonlimsAH = np.array([-122.72, -122.62])
latlimsAH = np.array([48.12, 48.18])

# Functionality copied from https://github.com/clawpack/geoclaw/blob/master/src/python/geoclaw/topotools.py#L873
land_cmap = plt.get_cmap('Greens_r')
sea_cmap = plt.get_cmap('Blues_r')
cmapPS = colormaps.add_colormaps((land_cmap, sea_cmap), data_limits=[-375,2500], data_break=0.0)
cmapAI = 'Blues_r'
cmapAH = 'Blues_r'

# levels to plot
levsPS = np.concatenate((np.arange(-375, 0, 25), np.arange(0,3000,500)))
levsAI = np.arange(-200, 20, 20)
levsAH = np.arange(-120, 15, 15)

# use basemap
basemapPS = Basemap(llcrnrlon=lonlimsPS[0], llcrnrlat=latlimsPS[0], 
                urcrnrlon=lonlimsPS[1], urcrnrlat=latlimsPS[1], 
                lat_0=latlimsPS.mean(), lon_0=lonlimsPS.mean(),
                projection='lcc', resolution='f',
                area_thresh=0.)
xPS, yPS = basemapPS(mat['lon_topo'], mat['lat_topo'])
xlimsAI, ylimsAI = basemapPS(lonlimsAI, latlimsAI)
xlimsAH, ylimsAH = basemapPS(lonlimsAH, latlimsAH)

# Make Puget Sound plot
fig = plt.figure(figsize=(16,16))
axPS = fig.add_subplot(111)
basemapPS.drawcoastlines(ax=axPS)
mappablePS = axPS.contourf(xPS, yPS, mat['z_topo'], cmap=cmapPS, levels=levsPS, zorder=2)
locator = MaxNLocator(6) # if you want no more than 10 contours
locator.create_dummy_axis()
locator.set_bounds(lonlimsPS[0], lonlimsPS[1])
pars = locator()
locator = MaxNLocator(6) # if you want no more than 10 contours
locator.create_dummy_axis()
locator.set_bounds(latlimsPS[0], latlimsPS[1])
mers = locator()
basemapPS.drawparallels(mers, dashes=(1, 1), linewidth=0.15, labels=[1,0,0,0], ax=axPS)#, zorder=3)
basemapPS.drawmeridians(pars, dashes=(1, 1), linewidth=0.15, labels=[0,0,0,1], ax=axPS)#, zorder=3)
cbPS = fig.colorbar(mappablePS, pad=0.015, aspect=35)
cbPS.set_label('Height/depth [m]')
# Label
axPS.text(0.8, 0.025, 'Puget Sound', transform=axPS.transAxes, color='0.15')

# Inset magnified plot of Admiralty Inlet
axAI = zoomed_inset_axes(axPS, 2, loc=1)
basemapPS.drawcoastlines(ax=axAI)
basemapPS.fillcontinents('darkgreen', ax=axAI)
mappableAI = axAI.contourf(xPS, yPS, mat['z_topo'], cmap=cmapAI, levels=levsAI)
axAI.set_xlim(xlimsAI)
axAI.set_ylim(ylimsAI)
# Inlaid colorbar
caxAI = fig.add_axes([0.582, 0.665, 0.011, 0.1])
cbAI = plt.colorbar(mappableAI, cax=caxAI, orientation='vertical')
cbAI.ax.tick_params(labelsize=12)
# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(axPS, axAI, loc1=2, loc2=4, fc="none", ec="0.3", lw=1.5, zorder=5)
# Label
axAI.text(0.41, 0.83, 'Admiralty\n      Inlet', transform=axAI.transAxes, color='0.15', fontsize=16)

# Inset magnified plot of Admiralty Head
axAH = zoomed_inset_axes(axPS, 9, loc=3)
basemapPS.drawcoastlines(ax=axAH)
basemapPS.fillcontinents('darkgreen', ax=axAH)
mappableAH = axAH.contourf(xPS, yPS, mat['z_topo'], cmap=cmapAH, levels=levsAH)
axAH.set_xlim(xlimsAH)
axAH.set_ylim(ylimsAH)
# Inlaid colorbar
caxAH = fig.add_axes([0.399, 0.116, 0.012, 0.15])
cbAH = plt.colorbar(mappableAH, cax=caxAH, orientation='vertical')
cbAH.ax.tick_params(labelsize=12)
# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(axPS, axAH, loc1=2, loc2=4, fc="none", ec="0.3", lw=1.5, zorder=5)
# Label
axAH.text(0.47, 0.92, 'Admiralty Head', transform=axAH.transAxes, color='0.15', fontsize=16)

# Save figure
fig.savefig('figures/domains.png', bbox_inches='tight')
# fig.show()


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