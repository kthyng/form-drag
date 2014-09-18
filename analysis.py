'''
Functions to do form drag analysis.
'''

import numpy as np
from matplotlib import delaunay
import netCDF4 as netCDF
import glob.glob as glob
import tracpy
from scipy import ndimage
import octant

# def interp2pt(xvar, yvar, zvar, var, xint, yint, zint):
#     '''
#     Interpolate a variable on a 3D field to input x,y,z locations.

#     Inputs:
#      xvar, yvar     The x, y locations of the variable var. 
#                     Should be size [NxM] where N is the y 
#                     dimension and M is the x dimension. Can be 
#                     unstructured.
#      zvar           The z locations of the variable var. 
#                     Should be size [FILL IN]. x and y can be 
#                     unstructured.
#      var            3D field to be interpolated. Should be size
#                     [FILL IN]
#      xint, yint     The x, y, z locations at which we want to 
#                     know var. Size [FILL IN]
#      zint           The x, y, z locations at which we want to 
#                     know var. Size [FILL IN]

#     Outputs:
#      varint             Interpolated variable
#     '''

#     ## Interpolate in x and y simulataneously, onto xint, yint, ##
#     ## carrying along the z dimension.                          ##

#     # Solve for x-y grid triangulation
#     tri = delaunay.Triangulation(xvar.flat, yvar.flat)

#     # loop through 

#     # Set up function for var
#     f = tri.nn_interpolator(var)

#     # Solve for var in x and y
#     varint = f(xint, yint)


#     # Interpolate 1D in z onto zint.


def run():
    '''
    Controlling script.
    '''

    # What files to use
    files = glob('ocean_his_000?.nc') # PLACE HOLDER FOR NOW

    # What time indices to use
    tinds = np.arange(200) # PLACE HOLDER FOR NOW

    # File objects
    nc = netCDF.MFDataset(files)

    # grid information
    loc = grid.nc # PLACEHOLDER
    grid = tracpy.inout.readgrid(loc)
    # xr = nc.variables['x_rho'][:] # [y,x], rho grid
    # yr = nc.variables['y_rho'][:] # [y,x], rho grid
    # dx = xr[0,1] - xr[0,0] # uniform grid spacing
    # dy = yr[1,0] - yr[0,0] # uniform grid spacing

    # Calculate required fields at native points on grid. Need:
    # rho, u, v, dh/dx, dh/dy, where h is the bathymetry
    rho = nc.variables['rho'][tinds,:,:,:] # [t,z,y,x], rho grid
    u = nc.variables['u'][tinds,:,:,:] # [t,z,y,x], u grid
    v = nc.variables['v'][tinds,:,:,:] # [t,z,y,x], v grid
    zeta = nc.variables['zeta'][tinds,:,:] # [t,y,x], ssh, rho grid
    h = nc.variables['h'][:,:] # [y,x], rho grid
    dhdy, dhdx = np.gradient(h, dx, dy) # [y,x], rho grid

    # Need depths in time
    zw = octant.depths.get_zw(grid['Vtransform'], grid['Vstretching'], 
            grid['km']+1, grid['theta_s'], grid['theta_b'], 
            h.T.copy(order='c'), grid['hc'], zeta=zeta, Hscale=3)
    zr = octant.depths.get_zrho(grid['Vtransform'], grid['Vstretching'], 
            grid['km'], grid['theta_s'], grid['theta_b'], 
            h.T.copy(order='c'), grid['hc'], zeta=zeta, Hscale=3)

    # Transect x,y locations in real space
    # NEED GRID TO SEE WHAT THESE ARE, AND NEED IN X,Y NOT LON/LAT
    # xt, yt

    # Convert transect locations from x,y space to grid index space
    xti, yti, dt = tracpy.tools.interpolate2d(xt, yt, grid, 'd_xy2ij')
    zti = 0.5 # PLACEHOLDER

    # Interpolate in 3D grid space to transect. The z direction can
    # be ignored for now because it will just go to the same part of 
    # the layer as from the original locations.
    rhot = ndimage.map_coordinates(rho, np.array([xi.flatten(), 
                                yi.flatten(), zi.flatten()]), 
                                order=3, mode='nearest',cval=0.).reshape(zi.shape)
    # SAME FOR U AND V
    # HOW TO GET DH/DX', THE ALONG-TRACK GRADIENT? PROJECT ONTO TRACK DIRECTION?

    # Calculate depths along transect (HOW TO DEAL WITH TIME HERE?)
    zwt = ndimage.map_coordinates(zw, np.array([xi.flatten(), 
                                yi.flatten(), zi.flatten()]), 
                                order=3, mode='nearest',cval=0.).reshape(zi.shape)
    zrt = ndimage.map_coordinates(zr, np.array([xi.flatten(), 
                                yi.flatten(), zi.flatten()]), 
                                order=3, mode='nearest',cval=0.).reshape(zi.shape)


    # Perform integrations, including critical depth z0
    # Along track

    # Vertically. Integrate below and above critical depth easily, and include 
    # a weighting of the percent of depth used for layers that cross z0


if __name__ == "__main__":
    run()    
