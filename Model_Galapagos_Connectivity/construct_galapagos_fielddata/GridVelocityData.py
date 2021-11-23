# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 10:03:16 2021

@author: quint
"""





import xarray as xr
import numpy as np

################################################### FUNCTIONS ######################################################
def getclosest_ij(lats,lons,latpt,lonpt):    
    """Function to find the index of the closest point to a certain lon/lat value."""
    dist_lat = (lats-latpt)**2          # find squared distance of every point on grid
    dist_lon = (lons-lonpt)**2
    minindex_lat = dist_lat.argmin()    # 1D index of minimum dist_sq element
    minindex_lon = dist_lon.argmin()
    return minindex_lat, minindex_lon   # Get 2D index for latvals and lonvals arrays from 1D index


#################################################### GRID AND VELOCITY DATA #########################################

PathToVelocity = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input\RGEMS_2010_Surf.nc'
PathToGrid = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input\RGEMS3_Surf_grid.nc'

#################################################### GET INDDICES OF DOMAIN OF INTEREST #########################################




griddata = xr.open_dataset(PathToGrid)
lonCentre = griddata['XC'].data
latCentre = griddata['YC'].data

#lonCorner, latCorner is the left lower corner of a grid
lonCorner = griddata['XG'].data  
latCorner = griddata['YG'].data

diff = latCorner - np.roll(latCorner,1)

lonCentre_grid = np.zeros((len(latCentre),len(lonCentre)))
latCentre_grid = np.zeros((len(latCentre), len(lonCentre)))

for row in range(len(lonCentre_grid[:,0])):
    lonCentre_grid[row,:] = lonCentre

for column in range(len(latCentre_grid[0,:])):
    latCentre_grid[:, column] = latCentre

lonCorner_grid = np.zeros((len(latCorner),len(lonCorner)))
latCorner_grid = np.zeros((len(latCorner), len(lonCorner)))

for row in range(len(lonCorner_grid[:,0])):
    lonCorner_grid[row,:] = lonCorner

for column in range(len(latCorner_grid[0,:])):
    latCorner_grid[:, column] = latCorner


domain = [-92, -88, -2, 2]

iy_min, ix_min = getclosest_ij(latCentre, lonCentre, domain[2], domain[0])
iy_max, ix_max = getclosest_ij(latCentre, lonCentre, domain[3], domain[1])

#################################################### COMPUTE VELOCITY AND GRID FIELDS/DATA #########################################

latCentredy = griddata['YC'][iy_min:iy_max].data
dx_g = (lonCentre[-1]-lonCentre[0])/len(lonCentre)
dy_g = (latCentredy[-1]-latCentredy[0])/len(latCentredy)

velocity_data = xr.open_dataset(PathToVelocity)
uvel = velocity_data['UVEL'][:,iy_min:iy_max,ix_min:ix_max].data
vvel = velocity_data['VVEL'][:,iy_min:iy_max,ix_min:ix_max].data

x_corners = lonCorner[ix_min:ix_max]
y_corners = latCorner[iy_min:iy_max]

x_corners_grid = np.zeros((len(y_corners),len(x_corners)))
y_corners_grid = np.zeros((len(y_corners), len(x_corners)))

for row in range(len(x_corners_grid[:,0])):
    x_corners_grid[row,:] = x_corners

for column in range(len(x_corners_grid[0,:])):
    y_corners_grid[:,column] = y_corners


#################################################### SAVE DATA #########################################

griddata.close() 
velocity_data.close()

np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\dx', dx_g)
np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\dy', dy_g)

np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\x_corners_grid', x_corners_grid)
np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\y_corners_grid', y_corners_grid)

np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\uvel', uvel)
np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\vvel', vvel)

np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\lonCentreGrid', lonCentre_grid)
np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\latCentreGrid', latCentre_grid)
np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\lonCornerGrid', lonCorner_grid)
np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\latCornerGrid', latCorner_grid)

np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\lonCentre', lonCentre)
np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\latCentre', latCentre)
np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\lonCorner', lonCorner)
np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\latCorner', latCorner)























