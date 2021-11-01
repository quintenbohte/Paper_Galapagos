# -*- coding: utf-8 -*-
"""
Created on Sat Oct 30 13:50:51 2021

@author: quint
"""

import numpy as np
import pandas as pd
from pathlib import Path
import os
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
#########################################DEFINE RELATIVE PATH#######################################

parent_dir = str(Path(os.path.abspath(__file__)).parents[2])

########################################## LOAD DATA ######################################

savename = 'long_run'
transitionmatrix = np.load(parent_dir + '\data\output\Transition_matrices\Tm_' +savename + '.npy')
gridsdataframe = pd.read_csv(parent_dir + '\data\input\galapagos_field_data\GridsDataFrame.csv')

########################################## DETERMINE TAKEOFF ###########################################


abs_numb_to_other_island = np.zeros((len(gridsdataframe)))

for release_grid in range(len(transitionmatrix[:,1])):
    release_island = gridsdataframe['island_number'][release_grid]
#    print('release island is:', release_island)
    for beach_grid in range(len(transitionmatrix[1,:])):
        
        beach_island = gridsdataframe['island_number'][beach_grid]
        print('beach island is:', beach_island)

        if release_island != beach_island:
            print('test')
            
            number_particles = transitionmatrix[release_grid, beach_grid]
            abs_numb_to_other_island[release_grid] += number_particles
        
perc_to_other_island = (abs_numb_to_other_island/1450)*100
        


############################################ PLOT TAKEOFF #################################################

scaling = 0.3
NetworkDataFrame = pd.DataFrame({'GridNumber':gridsdataframe['Grid_Number'], 'longitude':gridsdataframe['min_lon'], 'latitude':gridsdataframe['min_lat'],
                                 'abs_number_takeoff': abs_numb_to_other_island, 'perc_number_takeoff':perc_to_other_island})

plt.figure(figsize = (10,10))
ax = plt.axes(projection = ccrs.PlateCarree())
ax.set_extent([-92, -89, -1.5, 0.75], ccrs.PlateCarree())         
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.COASTLINE)
grid_lines = ax.gridlines(draw_labels=True, linestyle = '-.', color = 'black')
grid_lines.xformatter = LONGITUDE_FORMATTER
grid_lines.yformatter = LATITUDE_FORMATTER
plt.scatter(NetworkDataFrame['longitude'], NetworkDataFrame['latitude'], c = NetworkDataFrame['perc_number_takeoff'], s = NetworkDataFrame['abs_number_takeoff']*scaling)
plt.set_cmap('Reds')
plt.colorbar(orientation="horizontal").set_label('Betweenness Centrality (percentage most likely paths crossed)', fontsize = 15)



























