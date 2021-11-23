# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 11:34:25 2021

@author: quint
"""


###############################################################################
# Importing Packages                                                          #
###############################################################################

from datetime import timedelta as delta
from os import path
from glob import glob
import numpy as np
import pandas as pd
import dask
import math
import xarray as xr
from netCDF4 import Dataset
import warnings
import matplotlib.pyplot as plt
from operator import attrgetter

warnings.simplefilter('ignore', category=xr.SerializationWarning)

from parcels import AdvectionRK4
from parcels import Field
from parcels import FieldSet
from parcels import JITParticle
from parcels import ParticleFile
from parcels import ParticleSet
from parcels import Variable
from parcels import plotTrajectoriesFile
from parcels import rng as random


###############################################################################
# Used Functions                                                         #
###############################################################################
def getclosest_ij(lats,lons,latpt,lonpt):    
    """Function to find the index of the closest point to a certain lon/lat value."""
    dist_lat = (lats-latpt)**2          # find squared distance of every point on grid
    dist_lon = (lons-lonpt)**2
    minindex_lat = dist_lat.argmin()    # 1D index of minimum dist_sq element
    minindex_lon = dist_lon.argmin()
    return minindex_lat, minindex_lon   # Get 2D index for latvals and lonvals arrays from 1D index

def get_lon_lat_indices(lons,lats,lonloc,latloc):
    dist_lat = abs(lats - latloc)[:,0]
    dist_lon = abs(lons - lonloc)[0,:]

    min_y_index = np.argmin(dist_lat)
    min_x_index = np.argmin(dist_lon)
    
    return min_y_index, min_x_index


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
    GridSet = Dataset(grids)
    lon = GridSet.variables['XC'][:]
    lat = GridSet.variables['YG'][:]
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
    


def GetCoastGrids(LandMask):
    
    """
    This function will find all the CoastGrids, the first sea grids beside the land. 
    Input:
        - A landmask. 
    Output:
        - A 2D array with coastgrids = 1 and rest = 0.
    
    """
    
    """
    Define a coastline map. This map will be set to 1 on all coast cells.
    """
    
    
    CoastlineMap = np.zeros((LandMask.shape[0], LandMask.shape[1]))
    
    """
    We will use a nested loop to loop through all cells of the Landmask cell. What this loop basically does is,
    when a cell has a value of 1 (land), it will make all surrounding cells 1, so we create kind of an extra line of 
    grids around the landmask. In the end we will substract the  landmask from the mask which is created by the nested loop, 
    which result in only a mask with the coast grids. Notice, that when we're in the corner, upper, side, or lower row, and we
    meet a land cell, we should not make all surrounding cells 1. For example, we the lower left corner is a land grid, you should only make the inner cells 1. 
    """
    
    for i in range(LandMask.shape[0]-1):
        for j in range(LandMask.shape[1]-1):
            

            """
            We have nine if statements, four for the corners, four for the sides and one for the middle
            of the landmask. 
            """

            if i == 0 and j == 0: #upper left corner
                
                if LandMask[i,j] == 1:
                    
                    CoastlineMap[i,j] = 1
                    CoastlineMap[i,j+1] = 1
                    
                    CoastlineMap[i+1,j] = 1         
                    CoastlineMap[i+1, j+1] = 1
                    
                                
            elif i == 0 and j != 0 and j != LandMask.shape[1]-1: #upper row
                
                if LandMask[i,j] == 1:
                    
                    CoastlineMap[i,j] = 1
                    CoastlineMap[i,j-1] = 1
                    CoastlineMap[i,j+1] = 1
                    
                    CoastlineMap[i+1, j] = 1
                    CoastlineMap[i+1,j-1] = 1
                    CoastlineMap[i+1,j+1] = 1
             
                
            elif i == 0 and j == LandMask.shape[1]-1: #upper right corner
                
                if LandMask[i,j] == 1:
                
                    CoastlineMap[i,j] = 1
                    CoastlineMap[i,j-1] = 1
                    
                    CoastlineMap[i+1,j] = 1         
                    CoastlineMap[i+1, j-1] = 1
                    
            elif i != 0 and i != LandMask.shape[0]-1 and j == LandMask.shape[1]-1: #right row
                
                if LandMask[i,j] == 1:
                
                    CoastlineMap[i,j] = 1
                    CoastlineMap[i+1,j] = 1
                    CoastlineMap[i-1,j] = 1
                    
                    CoastlineMap[i, j-1] = 1
                    CoastlineMap[i+1,j-1] = 1
                    CoastlineMap[i-1,j-1] = 1
                    
            elif i == LandMask.shape[0]-1 and j == LandMask.shape[1]-1: #lower right corner
                
                if LandMask[i,j] == 1:
                
                    CoastlineMap[i,j] = 1
                    CoastlineMap[i,j-1] = 1
                    
                    CoastlineMap[i-1,j] = 1         
                    CoastlineMap[i-1, j-1] = 1
                    
            elif i == LandMask.shape[0]-1 and j != 0 and  j != LandMask.shape[1]-1: #lower row
                
                if LandMask[i,j] == 1:
                
                    CoastlineMap[i,j] = 1
                    CoastlineMap[i,j-1] = 1
                    CoastlineMap[i,j+1] = 1
                    
                    CoastlineMap[i-1, j] = 1
                    CoastlineMap[i-1,j-1] = 1
                    CoastlineMap[i-1,j+1] = 1
                
            
            elif i == LandMask.shape[0]-1 and j == 0: #lower left corner
                
                if LandMask[i,j] == 1:
                
                    CoastlineMap[i,j] = 1
                    CoastlineMap[i,j+1] = 1
                    
                    CoastlineMap[i+1,j] = 1         
                    CoastlineMap[i+1, j+1] = 1
            
            elif i != 0 and i != LandMask.shape[0]-1 and j == 0: #left row
                
                if LandMask[i,j] == 1:
                
                    CoastlineMap[i,j] = 1
                    CoastlineMap[i+1,j] = 1
                    CoastlineMap[i-1,j] = 1
                    
                    CoastlineMap[i, j+1] = 1
                    CoastlineMap[i+1,j+1] = 1
                    CoastlineMap[i-1,j+1] = 1
                    
            else:
            
                if LandMask[i,j] == 1:
                
                    CoastlineMap[i,j] = 1 #middle
                    CoastlineMap[i+1,j] = 1#lowermiddle
                    CoastlineMap[i-1,j] = 1#uppermiddle
                
                    CoastlineMap[i+1, j-1] = 1
                    CoastlineMap[i-1, j-1] = 1
                    CoastlineMap[i, j-1] =1
                
                    CoastlineMap[i+1, j+1] = 1
                    CoastlineMap[i-1, j+1] = 1
                    CoastlineMap[i, j+1] = 1
        
    
    
    """
    Here we substract the landmaks from the coastline mask, resulting in only
    the coastline. 
    """
    
    
    Coastgrids = CoastlineMap - LandMask
    
    return Coastgrids, CoastlineMap




def FindLocationCornerCoastGrid(landmask, coastgrids, fieldset, A = False):
    
    """
    This function creates a list with the longitudes and latitudes of lower left corner of the coast grids
    Notice the A argument. This is due to the inconvenient shape of this island. Will
    be explained in more detail later. 
    Input:
        - Landmask, coastgrid, fieldset
    Output:
        - List with longitudes and latitudes of the release locations of the particles. 
    """
    
    fU = fieldset.U
    
    
    """
    If we're using this function for island A, we set A = True. Then the longitude matrix will be
    shorten one row. This is needed because the island coast grids are determined in four seperate parts.
    However, the following lines of code will determine the longitudes and latitudes of the given fieldset. These will
    be used to set the correct longitudes and latitudes to coastgrids.
    """
    
    
    if A == False:
        lon = np.array(fU.lon[:]) 
        lat = np.array(fU.lat[:])
    
    if A == True:
        lon = np.array(fU.lon[:]) 
        lat = np.array(fU.lat[1:])
    
    
    
    """
    Here we make the lon and lat matrix, which will be filled with the longitude and latitudes 
    defined above. These lon lat matrices will have the same shape as the coastgrids mask.
    """
    lon_matrix = np.zeros((landmask.shape[0], landmask.shape[1]))
    lat_matrix = np.zeros((landmask.shape[0], landmask.shape[1]))
    
    for i in range(landmask.shape[0]):
        lon_matrix[i,:] = lon
    
    for j in range(landmask.shape[1]):
        lat_matrix[:,j] = lat
    
    """
    Finally, we will fill to list with the start longitude and latitudes. We're looping through
    the coastgrids mask and lon lat matrices, when a cell in the coastgrids mask is 1, we will append the corresponding lon and 
    lat (from the lon lat matrix we're also looping throug) in the list. 
    """
    startlon_list = []
    startlat_list = []
    
    for i in range(landmask.shape[0]):
        for j in range(landmask.shape[1]):
            if coastgrids[i,j] == 1:
                startlon_list.append(lon_matrix[i,j])
                startlat_list.append(lat_matrix[i,j])
                

    return startlon_list, startlat_list



def CreateDfWithFullGridInfo(releaselon, releaselat, IslandLetter, IslandNumber):
    
    """
    This function transforms the lists with release locations into a dataframe which contains:
        - min lon
        - max lon
        - min lat
        - max lat
        - Island Letter
        - Island number
        
    Input:
        - ReleaseLonList
        - ReleaseLatList
        - Island 
    """
    
    
    """
    Load dx and dy, to define grid sides from corner point
    """
    dx_g = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\dx.npy')
    dy_g = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\dy.npy')
    
    """
    Create dataframe of the two lists. We also sort the grids clockwise for every island. With the begin as the 
    most northern point of the island
    """
    ReleaseLocationsDf = pd.DataFrame({'lon': np.asarray(releaselon) - 0.5*dx_g, 'lat': np.asarray(releaselat) + 0.5*dy_g})
#    ReleaseLocationsDfCW = SortCoastGridsClockwise(ReleaseLocationsDf)
    ReleaseLocationsDf.index = np.arange(0,len(ReleaseLocationsDf))

    """
    Add the max lon, max lat, island number and island letter to the dataframe
    """
    ReleaseLocationsDf['max_lon'] = ReleaseLocationsDf['lon'].apply(lambda x: x + dx_g)
    ReleaseLocationsDf['max_lat'] = ReleaseLocationsDf['lat'].apply(lambda x: x - dy_g)
    ReleaseLocationsDf['Island_Letter'] = pd.DataFrame([IslandLetter]*len(ReleaseLocationsDf))
    ReleaseLocationsDf['Island_Number'] = pd.DataFrame([IslandNumber]*len(ReleaseLocationsDf))
    ReleaseLocationsDf['GridNumber'] = np.arange(1,len(ReleaseLocationsDf)+1)
    ReleaseLocationsDf.columns = ['min_lon', 'max_lat', 'max_lon', 'min_lat', 'Island_Letter', 'Island_Number', 'Grid_Number']
    ReleaseLocationsDf = ReleaseLocationsDf[['min_lon', 'min_lat', 'max_lon', 'max_lat', 'Island_Letter', 'Island_Number', 'Grid_Number']]

    return ReleaseLocationsDf

def CreateReleasePoints(points_on_longitude, points_on_latitude, grids):

    """
    This function can be used to determine the amount of particles we wanna release from the coast grids.
    input: 
        - point on the x-axis (points on longitude)
        - point on y-axis (points on latitude)
        - DataFrame with grid information (max lon, min lon, amax lat, min lat)
    
    output:
        - List with release longitudes and latitudes 
    """
    
    ReleasePointsLon = []
    ReleasePointsLat = []
    
    GridsCW_array = np.asarray(grids[['min_lon', 'min_lat', 'max_lon', 'max_lat']])
    
    for i in range(len(GridsCW_array)):
    
        lon_space = np.linspace(GridsCW_array[i,0], GridsCW_array[i,2], num = points_on_longitude+2 )
        lat_space = np.linspace(GridsCW_array[i,1], GridsCW_array[i,3], num = points_on_latitude+2 )
        
        
        lon_space_cor = lon_space[1:-1]
        lat_space_cor = lat_space[1:-1]
        
        for j in lon_space_cor:
            for k in lat_space_cor:
                                    
                ReleasePointsLon.append(j)
                ReleasePointsLat.append(k)
            
    return ReleasePointsLon, ReleasePointsLat
###############################################################################
# Define Release Locations of Particles for every Island                      #
###############################################################################


#Loading dx, dy, grid and velocity data        
dx_g = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\dx.npy')
dy_g = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\dy.npy')
PathToVelocity = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input\RGEMS_2010_Surf.nc'
PathToGrid = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input\RGEMS3_Surf_grid.nc'


domain = [-92, -88, -2, 2]

#create fieldset for island A1
fieldset = fieldset(domain, PathToVelocity, PathToGrid)


#create landmask for island A1
Landmask_galapagos = CreateLandmask(fieldset)
Landmask_galapagos = Landmask_galapagos[1:, :]
Landmask_galapagos[88,:] = 0
Landmask_galapagos[:,119] = 0

#Create coast grid mask
CoastGrids, Coastline = GetCoastGrids(Landmask_galapagos)

#find longitudes and latitudes of the lower left corner of the coastgrids
releaselon, releaselat = FindLocationCornerCoastGrid(Landmask_galapagos, CoastGrids, fieldset, A=True)
releaselon = releaselon + 0.5*dx_g
releaselat = releaselat + 0.5*dy_g


GridsDataFrame = CreateDfWithFullGridInfo(releaselon, releaselat, IslandLetter = 'X', IslandNumber = 1)

np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\ReleaseLocations\ReleaseLon', releaselon)
np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\ReleaseLocations\ReleaseLat', releaselat)
np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\Coastgrids', CoastGrids )
np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\Coastline', Coastline )
GridsDataFrame.to_csv(r'C:\Users\quint\Documents\Quinten_studie\Master_Thesis\Data\Output_Data\ReleaseLocations\GridsDataFrame.csv', index = False)






































