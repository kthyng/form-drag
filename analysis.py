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
import bisect
import plotting


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

    # save transect lat/lon for future use
    np.savez('transect_lonlat.npz', lont=lont, latt=latt)

    return lont, latt, theta, dd


def get_vars(nc, tinds, grid):
    '''
    Get and calculate some necessary fields for calculations.

    Inputs:
     nc         netCDF file object
     tinds      Time indices
     grid       Grid dictionary from TracPy

    Outputs:
     rho        Density with 1000 added in to be about 1023 [kg/m^3]
     u,v        x and y direction velocities
     dhdx       x-component of bathymetric gradient - only need one to 
                get along-transect component
     zr         Depths at vertical rho-grid points
     dzr        Layer thicknesses of the vertical layers
    '''

    ## Read in variables ##
    rho = nc.variables['rho'][tinds,:,:,:] # [t,z,y,x], rho grid
    u = nc.variables['u'][tinds,:,:,:] # [t,z,y,x], u grid
    v = nc.variables['v'][tinds,:,:,:] # [t,z,y,x], v grid
    h = nc.variables['h'][:,:] # [y,x], rho grid
    dx = 1./grid['pm'][0,0]
    dy = 1./grid['pn'][0,0]
    dhdy, dhdx = np.gradient(h, dx, dy) # [y,x], rho grid


    ## Calculate depths in time ##

    # Need ssh for calculating the depths
    zeta = nc.variables['zeta'][tinds,:,:] # [t,y,x], ssh, rho grid

    zw = np.empty((tinds.size, grid['km']+1, h.shape[0], h.shape[1]))
    zr = np.empty((tinds.size, grid['km'], h.shape[0], h.shape[1]))
    # loop through times
    for i,tind in enumerate(tinds):
        zw[i,:] = octant.depths.get_zw(grid['Vtransform'], grid['Vstretching'], 
                grid['km']+1, grid['theta_s'], grid['theta_b'], 
                h, grid['hc'], zeta=zeta[tind,:,:], Hscale=3)
        zr[i,:] = octant.depths.get_zrho(grid['Vtransform'], grid['Vstretching'], 
                grid['km'], grid['theta_s'], grid['theta_b'], 
                h, grid['hc'], zeta=zeta[tind,:,:], Hscale=3)

    # Get layer thicknesses
    dzr = zw[:,1:,:,:] - zw[:,:-1,:,:]

    return rho+1000., u, v, dhdx, zr, dzr


def interp2transect(sizes, rho, u, v, zr, dzr, dhdx, h, pts, theta):
    '''
    Interpolate from planar space to along the transect (x' direction) 
    in x,y, often bringing the time and vertical dimensions along 
    (which are not interpolated in).

    Inputs:
     sizes          [time,z,x'] array sizes
     rho, u, v, zr,
     dzr, dhdx, h   Fields to be interpolated
     pts            Points in grid space to find the fields at
     theta          Angle of the transect line

    Outputs:
     fields at the transect locations
    '''

    # Interpolate in 3D grid space to transect. The z direction can
    # be ignored for now because it will just go to the same part of 
    # the layer as from the original locations.
    # x' along-transect
    # [t,z(grid),y,x] --> [t,z(grid),x'] -->[t,z(real),x']
    # First, loop in time and z, then interpolate in x and y while leaving z
    rhot = np.empty(sizes)
    ut = np.empty(sizes)
    vt = np.empty(sizes)
    zrt = np.empty(sizes)
    dzrt = np.empty(sizes) # layer thicknesses
    # pdb.set_trace()
    for i in xrange(sizes[0]): # loop through time
        for k in xrange(sizes[1]): # loop through vertical layers
            # have to manipulate mask for this to work out
            rhotemp = rho[i,k,:,:].data
            ind = rhotemp>2000 # catch masked values
            rhotemp[ind] = 0.
            rhot[i,k,:] = ndimage.map_coordinates(rhotemp, pts, order=1, mode='nearest',cval=0.)
            utemp = u[i,k,:,:].data
            ind = utemp>100 # catch masked values
            utemp[ind] = 0.
            ut[i,k,:] = ndimage.map_coordinates(utemp, pts, order=1, mode='nearest',cval=0.)
            vtemp = v[i,k,:,:].data
            ind = vtemp>100 # catch masked values
            vtemp[ind] = 0.
            vt[i,k,:] = ndimage.map_coordinates(vtemp, pts, order=1, mode='nearest',cval=0.)
            # these gives the zr value along the transect, in depth
            zrt[i,k,:] = ndimage.map_coordinates(zr[i,k,:,:], pts, order=1, mode='nearest',cval=0.)
            dzrt[i,k,:] = ndimage.map_coordinates(dzr[i,k,:,:], pts, order=1, mode='nearest',cval=0.)
    # pdb.set_trace()
    # don't need to loop for bathymetry
    dhdxt = ndimage.map_coordinates(dhdx, pts, order=1, mode='nearest',cval=0.)
    # dhdyt = ndimage.map_coordinates(dhdy, pts, order=1, mode='nearest',cval=0.)

    # Calculate along-track direction (transect goes from left to right in the domain)
    dhdt = dhdxt/np.cos(theta)

    # depths along the transect
    ht = ndimage.map_coordinates(h, pts, order=1, mode='nearest',cval=0.)

    return rhot, ut, vt, zrt, dzrt, dhdt, ht


def form_drag(z0, dd, HT, sizes, zrt, dzrt, rhot, dhdt):
    '''
    Calculate the surface and internal form drag terms.

    Inputs:
     z0         Critical depth between internal and surface form drag definitions
     dd         Distance along transect between points [m]
     HT         Height of bump [m]
     sizes      [time,x'] array sizes
     zrt,dzrt   Depths/layer thicknesses along transect
     rhot       Densities along transect
     dhdt       Bathymetry gradient along transect

    Outputs:
     Dfi        Internal form drag in time [N/m^2]
     Dfs        Surface form drag in time [N/m^2]
    '''

    # Perform integrations:

    # --- Integrate vertically. Integrate below and above critical depth easily, and include 
    # a weighting of the percent of depth used for layers that cross z0
    # [t,z,x'] --> [t,x']
    Dfit = np.empty(sizes) # internal form drag
    Dfst = np.empty(sizes) # surface form drag
    # pdb.set_trace()
    for i in xrange(sizes[0]): # loop in time
        for j in xrange(sizes[1]): # loop along transect
            iz0 = bisect.bisect(zrt[i,:,j], z0) - 1 # index just below critical depth in depth
            # internal
            Dfit[i,j] = (rhot[i,:iz0,j]*dzrt[i,:iz0,j]).sum(axis=0) # up to the layer below the layer with z0
            dz = ((z0-zrt[i,iz0,j])/dzrt[i,iz0,j])*abs(z0) # fraction of layer below z0 * z0 to get in meters
            Dfit[i,j] = Dfit[i,j] + rhot[i,iz0,j]*dz # addition of part of layer that includes z0
            # surface
            Dfst[i,j] = (rhot[i,iz0+1:,j]*dzrt[i,iz0+1:,j]).sum(axis=0) # the layer above the layer with z0 and above
            dz = ((zrt[i,iz0+1,j]-z0)/dzrt[i,iz0,j])*abs(z0) # fraction of layer above z0 * z0 to get in meters
            Dfst[i,j] = Dfst[i,j] + rhot[i,iz0,j]*dz # addition of part of layer that includes z0

    # --- Integrate along track, [t,x'] --> [t]
    # Note that Dfi and Dfs have already been divided by W, the width of the bump
    g = 9.81 # gravity
    Dfi = -g * (Dfit*dhdt[np.newaxis,:].repeat(sizes[0],axis=0)).sum(axis=1)*dd
    Dfs = -g * (Dfst*dhdt[np.newaxis,:].repeat(sizes[0],axis=0)).sum(axis=1)*dd

    # Divide the form drags by the height of the bump, HT
    Dfi = Dfi/HT
    Dfs = Dfs/HT

    return Dfi, Dfs


def bottom_friction(nc, ut, vt, theta, dd, HT):
    '''
    Calculate the bottom friction.

    Inputs:
     nc         netCDF file object
     ut, vt     u,v velocities along transect [t,z,x']
     theta      Angle of transect line
     dd         Distance along transect between points [m]
     HT         Height of bump [m]

    Outputs:
     Dbbl       Drag due to bottom friction [N/m^2]
    '''

    rho0 = nc.variables['rho0'][:] # background density
    CD = nc.variables['rdrg2'][:] # is this correct???
    st = np.sqrt(ut[:,0,:]**2+vt[:,0,:]**2)
    taubxt = rho0 * CD * ut[:,0,:] / st # take velocity at bottom/seabed
    # taubyt = rho0 * CD * vt / s

    # Calculate along-track direction (transect goes from left to right in the domain)
    taubt = taubxt/np.cos(theta)

    # Then average along the track, multiply by W, and divide by HT
    # [t,x'] --> [t]
    Dbbl = taubt.mean(axis=1) # along-track average
    W = dd*ut.shape[2] # width of transect/bump
    Dbbl = Dbbl*W/HT

    return Dbbl


def plot_drag(tinds, Dfi, Dfs, Dbbl, ut):
    '''
    Plot the components of drag in time. Follow Edwards et al paper layout.
    '''

    fig = plt.figure(figsize=(12,10))

    # Plot form drag components
    ax1 = fig.add_subplot(3,1,1)
    ax1.plot(tinds, Dfi, 'k', lw=3)
    ax1.plot(tinds, Dfs, 'k', lw=2, alpha=0.6)
    ax1.set_ylabel('D_{FORM} (WH_T)^{-1} [Nm^{-2}]')
    ax1.legend(('Internal', 'Surface'))

    # Plot bottom friction component
    ax2 = fig.add_subplot(3,1,2)
    ax2.plot(tinds, Dbbl, 'k', lw=3)
    ax2.set_ylabel('D_{BBL} (WH_T)^{-1} [Nm^{-2}]')

    # Plot the velocity signal - start with U but change to PCA version soon
    # UPDATE THIS
    ax3 = fig.add_subplot(3,1,3)
    ax3.plot(tinds, ut[:,-1,0], 'k', lw=3)
    ax3.set_ylabel('u [ms^{-1}]')
    ax3.set_xlabel('Time')

    # ADD ALL THE LABELS

    fig.savefig('figures/drag.pdf', bbox_inches='tight')


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
    grid = tracpy.inout.readgrid(files[0], llcrnrlon=-122.8, llcrnrlat=47.9665, 
                                            urcrnrlon=-122.54, urcrnrlat=48.227, 
                                            lat_0=48, lon_0=-122.7, usebasemap=True, res='i')

    # Transect x,y locations in real space and angle of transect
    lont, latt, theta, dd = plot_domain(grid, nc)

    # Calculate required fields at native points on grid. 
    # Need: rho, u, v, dh/dx, dh/dy, where h is the bathymetry
    rho, u, v, dhdx, zr, dzr = get_vars(nc, tinds, grid)

    # Convert transect locations from x,y space to grid index space
    xgt, ygt, dt = tracpy.tools.interpolate2d(lont, latt, grid, 'd_ll2ij')
    zgt = np.arange(rho.shape[1]) # keep output at the layer locations in grid space

    # Array lengths: time, vertical layers, transect lengths
    tl = tinds.size; zl = zgt.size; xl = xgt.size

    # pdb.set_trace()

    # Find fields at transect locations
    # tracpy grid stuff is transposed due to tracmass
    rhot, ut, vt, zrt, dzrt, dhdt, ht = interp2transect((tl, zl, xl), rho, u, v, zr, dzr, 
                                                        dhdx, np.asarray(grid['h'].T, order='C'), np.array([ygt, xgt]), theta)

    # Calculate height of bump
    HT = ht.max()-ht.min()
    # # check the depths along the transect
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(np.arange(xgt.size)*dd, -ht, '0.2', lw=3)
    # ax.set_xlabel('Distance along transect [m]')
    # ax.set_ylabel('Depth [m]')
    # fig.show()

    # pdb.set_trace()

    # Add transect to larger domain plot and depths and save
    plotting.ps(dosave=True, fname='figures/domains-transect.png', lont=lont, latt=latt, ht=ht, dd=dd)

    # Integrate to get form drag terms
    z0 = -10 # critical depth for integrating surface vs. internal components
    Dfi, Dfs = form_drag(z0, dd, HT, (tl,xl), zrt, dzrt, rhot, dhdt)

    # Calculate the bottom friction from ROMS
    Dbbl = bottom_friction(nc, ut, vt, theta, dd, HT)

    # Plot drag results
    plot_drag(tinds, Dfi, Dfs, Dbbl, ut)


if __name__ == "__main__":
    run() 

