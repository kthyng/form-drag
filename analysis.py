'''
Functions to do form drag analysis.
'''

import numpy as np
from matplotlib import delaunay
import netCDF4 as netCDF
import glob
import tracpy
from scipy import ndimage
import octant
import pdb
import matplotlib.pyplot as plt
import tracpy.plotting
import tracpy.calcs


def plot_domain(grid, nc, dd=65):
    '''
    Plot the domain and then choose end points of a transect.

    Inputs:
     grid       TracPy grid dictionary
     nc         netCDF file object
     dd         spacing along transect (65 meters)

    Outputs:
     lont,latt  Longitude and latitude along the transect
     theta      Angle of the transect
     dd         Distance between points along the transect
    '''

    # plot magnification of Admiralty Head
    import plotting
    plotting.ai(grid, 'in') 

    # Have user select end points of transect
    print 'left click the two end points of the desired transect, end by closing plot.'
    print 'Choose points from upper left to lower right.'
    pts = plt.ginput(n=0, timeout=0) # pts in projected space
    x = []; y = []
    [x.append(pts[i][0]) for i in xrange(len(pts))]
    [y.append(pts[i][1]) for i in xrange(len(pts))]

    # x, y are the transect end points in projected space and lon/lat in geographic
    lon, lat = grid['basemap'](x, y, inverse=True)
    # endpoints in grid space
    xg, yg, _ = tracpy.tools.interpolate2d(x, y, grid, 'd_xy2ij')

    # Print depth at the end points since we want it to match
    depths = tracpy.calcs.Var(xg+1, yg+1, 0., 'h', nc) # to fortran indexing
    # MAKE THIS INTO AN ASSERTION
    ddiff = depths[1]-depths[0]
    if ddiff>2:
        print 'The end points are more than two meters apart in depth (' + str(ddiff) + '); try again'

    # Create line from end points
    m = (y[1] - y[0])/(x[1] - x[0]) # slope of transect
    # b = y[1] - m*x[1] # slope intercept
    theta = np.arctan(m) # angle of transect
    dxt = dd*np.cos(theta) # x increment, which has been projected from along the transect
    dyt = dd*np.sin(theta) # y increment, which has been projected from along the transect
    xt = np.arange(x[0], x[1], dxt)
    # yt = m*xt + b
    #or
    yt = np.arange(y[0], y[1], dyt)
    # pdb.set_trace()
    lont, latt = grid['basemap'](xt, yt, inverse=True)

    # Add transect to larger domain plot and save
    plotting.ps(dosave=True, fname='figures/domains-transect.png', lont=lont, latt=latt)

    return lont, latt, theta, dd



def run():
    '''
    Controlling script.
    '''

    # What files to use
    lochis = '/Volumes/Emmons/ai65/OUT/'
    files = glob.glob(lochis + 'ocean_his_00??.nc') # PLACE HOLDER FOR NOW

    # File objects
    nc = netCDF.MFDataset(files)

    # What time indices to use
    tinds = np.arange(9) # PLACE HOLDER FOR NOW

    # grid information
    # loc = grid.nc # PLACEHOLDER

    grid = tracpy.inout.readgrid(files[0], llcrnrlon=-122.8, llcrnrlat=47.9665, 
                                            urcrnrlon=-122.54, urcrnrlat=48.227, 
                                            lat_0=48, lon_0=-122.7, usebasemap=True, res='i')
    dx = 1./grid['pm'][0,0]
    dy = 1./grid['pn'][0,0]

    # Transect x,y locations in real space and angle of transect
    lont, latt, theta, dd = plot_domain(grid, nc)

    # Calculate required fields at native points on grid. Need:
    # rho, u, v, dh/dx, dh/dy, where h is the bathymetry
    rho = nc.variables['rho'][tinds,:,:,:] # [t,z,y,x], rho grid
    u = nc.variables['u'][tinds,:,:,:] # [t,z,y,x], u grid
    v = nc.variables['v'][tinds,:,:,:] # [t,z,y,x], v grid
    zeta = nc.variables['zeta'][tinds,:,:] # [t,y,x], ssh, rho grid
    h = nc.variables['h'][:,:] # [y,x], rho grid
    dhdy, dhdx = np.gradient(h, dx, dy) # [y,x], rho grid

    # Need depths in time
    # grid['theta_s'] = 0.001 # is actually 0
    # loop through times
    zw = np.empty((tinds.size, grid['km']+1, h.shape[0], h.shape[1]))
    zr = np.empty((tinds.size, grid['km'], h.shape[0], h.shape[1]))
    for i,tind in enumerate(tinds):
        zw[i,:] = octant.depths.get_zw(grid['Vtransform'], grid['Vstretching'], 
                grid['km']+1, grid['theta_s'], grid['theta_b'], 
                h, grid['hc'], zeta=zeta[tind,:,:], Hscale=3)
        zr[i,:] = octant.depths.get_zrho(grid['Vtransform'], grid['Vstretching'], 
                grid['km'], grid['theta_s'], grid['theta_b'], 
                h, grid['hc'], zeta=zeta[tind,:,:], Hscale=3)

    # Convert transect locations from x,y space to grid index space
    xgt, ygt, dt = tracpy.tools.interpolate2d(lont, latt, grid, 'd_ll2ij')
    zgt = np.arange(rho.shape[1]) # keep output at the layer locations in grid space

    # Interpolate in 3D grid space to transect. The z direction can
    # be ignored for now because it will just go to the same part of 
    # the layer as from the original locations.
    # x' along-transect
    # [t,z(grid),y,x] --> [t,z(grid),x'] -->[t,z(real),x']
    # First, loop in time and z, then interpolate in x and y while leaving z
    rhot = np.empty((tinds.size, zgt.size, xgt.size))
    ut = np.empty((tinds.size, zgt.size, xgt.size))
    vt = np.empty((tinds.size, zgt.size, xgt.size))
    zwt = np.empty((tinds.size, zgt.size, xgt.size))
    zrt = np.empty((tinds.size, zgt.size, xgt.size))
    for i,tind in enumerate(tinds): # loop through time
        for k in xrange(zgt.size): # loop through vertical layers
            # have to manipulate mask for this to work out
            rhotemp = rho[i,k,:,:].data
            ind = rhotemp>100 # catch masked values
            rhotemp[ind] = 0.
            rhot[i,k,:] = ndimage.map_coordinates(rhotemp, np.array([ygt, xgt]), 
                                    order=3, mode='nearest',cval=0.)
            utemp = u[i,k,:,:].data
            ind = utemp>100 # catch masked values
            utemp[ind] = 0.
            ut[i,k,:] = ndimage.map_coordinates(utemp, np.array([ygt, xgt]), 
                                    order=3, mode='nearest',cval=0.)
            vtemp = v[i,k,:,:].data
            ind = vtemp>100 # catch masked values
            vtemp[ind] = 0.
            vt[i,k,:] = ndimage.map_coordinates(vtemp, np.array([ygt, xgt]), 
                                    order=3, mode='nearest',cval=0.)
            # these gives the zw/zr value along the transect, in depth?
            zwt[i,k,:] = ndimage.map_coordinates(zw[i,k,:,:], np.array([ygt, xgt]), 
                                    order=3, mode='nearest',cval=0.)
            zrt[i,k,:] = ndimage.map_coordinates(zr[i,k,:,:], np.array([ygt, xgt]), 
                                    order=3, mode='nearest',cval=0.)

    # don't need to loop for bathymetry
    dhdxt = ndimage.map_coordinates(dhdx, np.array([ygt, xgt]), order=3, mode='nearest',cval=0.)
    dhdyt = ndimage.map_coordinates(dhdy, np.array([ygt, xgt]), order=3, mode='nearest',cval=0.)

    # Calculate along-track direction (transect goes from left to right in the domain)
    dhdt = dhdxt/np.cos(theta)

    # Perform integrations: TEST UNITS

    # --- Integrate along track
    g = 9.81
    Dalong = -g * (rhot*dhdt).sum(axis=1)*dd # common along-track component

    # --- Integrate vertically. Integrate below and above critical depth easily, and include 
    # a weighting of the percent of depth used for layers that cross z0
    pdb.set_trace()
    Dfi = Dalong # internal form drag
    Dfs = Dalong # surface form drag

if __name__ == "__main__":
    run()    
