'''
Plot bathymetry from Admiralty Inlet 65 meter simulation.
'''

import scipy.io
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

# download bathymetry, which can be found at: http://figshare.com/preview/_preview/1165560 (27.3MB)

# Read in bathymetry
mat = scipy.io.loadmat('cascadia_gridded.mat')

# x and y limits for this plot
lonlims = [-122.8, -122.55]
latlims = [47.9665, 48.227]
# lonlims = [-122.72, -122.65]
# latlims = [48.135, 48.18]

# # Grid for this plot
# lon = np.linspace(lonlims[0], lonlims[1], 1000)
# lat = np.linspace(latlims[0], latlims[1], 1000)
# [LON,LAT] = np.meshgrid(lon,lat)

# # Interpolate topo/bathy data to this grid
# Z = interpolate.interp2d(mat['lon_topo'], mat['lat_topo'], mat['z_topo'], LON, LAT)

# levels to plot
levs = np.linspace(-200, 200, 20)

# Make plot
fig = plt.figure()
ax = fig.add_subplot(111)
mappable = ax.contourf(mat['lon_topo'], mat['lat_topo'], mat['z_topo'], cmap='ocean_r', levels=levs)
ax.set_xlim(lonlims)
ax.set_ylim(latlims)
fig.colorbar(mappable)
# ax.pcolormesh(LON, LAT, Z, cmap='ocean')
# ax.plot(-122.6855,48.1515,'r.','markersize',20)
# % plot([-122.7128 -122.6699 -122.6699 -122.7128 -122.7128],...
# %     [48.1425 48.1425 48.1703 48.1703 48.1425],'k','linewidth',2)

# % Make plot pretty with 20 meter contour intervals
# demcmap('inc',Z,4) % from mapping toolbox

# % Add nice details to plot
# set(gcf,'position',[262   -13   738   697]) % Plot size
# axis equal
# axis tight
# ylabel(colorbar,'m','fontsize',20,'fontweight','bold') % Colorbar label
# set(gca,'fontsize',20,'fontweight','bold') % Label fonts

# % Save plot
# fname = 'bathy_ai_ahead_normal';
# saveas(gcf,fname,'fig')
# savefig(fname,'png','-r100')
# savefig(fname,'pdf','-r100')

