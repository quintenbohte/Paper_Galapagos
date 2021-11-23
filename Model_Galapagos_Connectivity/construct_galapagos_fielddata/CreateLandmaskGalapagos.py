# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 16:30:53 2021

@author: quint
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 16:45:14 2021

@author: quint
"""

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from parcels import FieldSet
import xarray as xr
import warnings
warnings.simplefilter('ignore', category=xr.SerializationWarning)
from glob import glob
import netCDF4 as nc
from netCDF4 import Dataset


###############################################################################
# Used Functions                                                              #
###############################################################################

def get_lon_lat_indices(lons,lats,lonloc,latloc):
    dist_lat = abs(lats - latloc)[:,0]
    dist_lon = abs(lons - lonloc)[0,:]

    min_y_index = np.argmin(dist_lat)
    min_x_index = np.argmin(dist_lon)
    
    return min_y_index, min_x_index


def getclosest_ij(lats,lons,latpt,lonpt):    
    """Function to find the index of the closest point to a certain lon/lat value."""
    dist_lat = (lats-latpt)**2          # find squared distance of every point on grid
    dist_lon = (lons-lonpt)**2
    minindex_lat = dist_lat.argmin()    # 1D index of minimum dist_sq element
    minindex_lon = dist_lon.argmin()
    return minindex_lat, minindex_lon   # Get 2D index for latvals and lonvals arrays from 1D index


def fieldset(domain, velocities, grids):
    
    from netCDF4 import Dataset
    from parcels import FieldSet
    """ Creating a Fieldset
    This function will create a fieldset from two netcdf files. One containing the velocities and the
    other one containing the grid information. As input, it needs the domain boundaries of the
    fieldset (Galapagos domain in my case), a path to the velocitiy netcdf file (velocities) and
    a path to the grid netcdf file (grids).
    """
    
    """
    The first set of lines will define the indices of the longitude and latitude used for the fieldset.
    So basically we will set the boundaries for the domain. It uses the get_lon_lat function, which is 
    defined above. 
    """
#    GridSet = Dataset(grids)
#    lon = GridSet.variables['XC'][:]
#    lat = GridSet.variables['YC'][:]
#    ymin, xmin = get_lon_lat_indices(lon, lat, domain[0], domain[2])
#    ymax, xmax = get_lon_lat_indices(lon, lat, domain[1], domain[3])
#    indices_MITgcm = {'lon': range(xmin,xmax), 'lat': range(ymin,ymax)}
#    
    
    data_in = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input'
    dfile = Dataset(data_in+'\RGEMS3_Surf_grid.nc')
    lon = dfile.variables['XC'][:]
    lat = dfile.variables['YC'][:]
    iy_min, ix_min = getclosest_ij(lat, lon, domain[2], domain[0])
    iy_max, ix_max = getclosest_ij(lat, lon, domain[3], domain[1])
    indices_MITgcm = {'lon': range(ix_min,ix_max), 
           'lat': range(iy_min,iy_max)}
        
        
    """
    The next lines of code will define the dictionaries needed for the fieldset. The first dictionary
    conatins the velocity and grid files. The second one contains the variable names and the thord one
    contains the dimension names. 
    """
    VelocityField = sorted(glob(velocities))
    GridField = glob(grids)
    files_MITgcm = {'U': {'lon': GridField, 'lat': GridField, 'data': VelocityField},                                             
                'V': {'lon': GridField, 'lat': GridField, 'data': VelocityField}}       
    variables_MITgcm = {'U': 'UVEL', 'V': 'VVEL'}
    dimensions_MITgcm = {'lon': 'XG', 'lat': 'YG', 'time': 'time'}
    
    
    """
    This last line of code creates the final fieldset
    """
    fieldset = FieldSet.from_mitgcm(files_MITgcm,
                                           variables_MITgcm, 
                                           dimensions_MITgcm,
                                           indices_MITgcm)
    
    return fieldset
    

def ChangeValues(Landmask, value_to_change, value_to_change_to):
    SideArray = np.zeros((len(Landmask[:,0]),len(Landmask[0,:])))
    for i in range(len(Landmask[:,0])):
        for j in range(len(Landmask[0,:])):
            
            value = Landmask[i,j]
            if value == value_to_change:
                SideArray[i,j] = value_to_change_to
                
    return SideArray
        




def CreateLandmask(Fieldset, test = False):
    
    """This function creates a Landmask of the fieldset    
    Land = 1 and Ocean = 0. 
    Input: fieldset, test =True/False
    Output: Landmask (land =1 , Ocean = 0)
    """
    
    
    """
    This first set of lines creates a numpy array with u velocities and a numpy
    array with v velocities. First we get the U and V fields from the dataset. Then
    we compute a time chunk, which is needed because of the dataset. Then we only
    take the first slice of the U and V field (we do not need more for finding the land
    and ocean grids). As last we make an empty array which will be filled with zeros and 
    ones.
    """
    fU = Fieldset.U
    fV = Fieldset.V
    Fieldset.computeTimeChunk(fU.grid.time[0], 1) 
    uvel_mask_c = fU.data[0,:,:] 
    vvel_mask_c = fV.data[0,:,:]
#    vvel_mask_c = np.roll(vvel_mask_c, 1, axis = 0)
    landmask =  np.zeros((uvel_mask_c.shape[0], uvel_mask_c.shape[1]))
    
    """
    The first loop checks the value of the u and v velocitites. Notice that we get the
    values of two adjacent grid, since we're working with a C-grid.
    Visualizations of velocities in the C-grids(see below). So for a grid to be flagged identified
    as a land grid two U velocities and 2 V velocities need to be zero. The first loop makes all
    ocean grids 1 and land grids 0.            
    ____ ____  ____  ____
    |    V    |     V     |           
    |         |           | 
    U    T    U     T     U
   |          |           |  
   |____V____|_____V_____|  
    """
    
    for i in range (len(landmask[:,0])-1):
        for j in range (len(landmask[0,:])-1):
            u1 = uvel_mask_c[i,j]

            u2 = uvel_mask_c[i,j+1]

            v1 = vvel_mask_c[i,j]

            v2 = vvel_mask_c[i+1,j]

            if u1 != 0 or u2 != 0 or v1 != 0 or v2 != 0:
                landmask[i,j] = 1
               
                
    """
    Change all zero to 1 and rest 0. since we want the land grids to be 1 and ocean
    grids to be 0. 
    """
            
    landmask = ChangeValues(landmask,0,1)        
    
    """
    The created landmask needs to be shifted upwards one grid. We will
    use the numpy roll function to do this.
    """
    
    if test == True:
        plt.figure()
        plt.imshow(landmask)
        plt.colorbar()
    
    return landmask
    

###############################################################################
# Creating Landmask                                                           #
###############################################################################


PathToVelocity = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Model_local_pc\data\input\RGEMS_surf\RGEMS_surf_all_years\*'
PathToGrid = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Model_local_pc\data\input\grid_data\RGEMS3_Surf_grid.nc'
domain = [-96, -84, -6, 6]

fieldset = fieldset(domain, PathToVelocity, PathToGrid)

Landmask_galapagos = CreateLandmask(fieldset)
#Landmask_galapagos[88,:] = 0
#Landmask_galapagos[:,119] = 0

plt.figure()
plt.imshow(Landmask_galapagos)
plt.colorbar()

#Save Landmask
np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\LandmaskLarge', Landmask_galapagos)






