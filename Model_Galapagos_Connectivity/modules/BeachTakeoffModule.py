# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 15:04:46 2021

@author: quint
"""


import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pathlib import Path
import os
parent_dir = str(Path(os.path.abspath(__file__)).parents[1])


#This class contains the metthods, which execute the beach and takeoff analysis. 

#TAKE OFF LOCATION:
#What is meant with take-off locations is: the number of particles that take-off from a grid that eventually end up on another island. So when a 
#particle beaches at the same island it is released from it is not taken into account in ths take-off metric for the node. 

#BEACH LOCATIONS:
#this is the opposite of the take off locations. Here compute the number of particles that beaches at a grid, but only take into account the 
#particle that are released from another island. 

class BeachTakeOffAnalysis():
    
    #Here we define the internal variables
    def __init__(self, TransitionMatrix):

        self.TransitionMatrix = TransitionMatrix
        self.GridDataFrame = pd.read_csv(parent_dir + '\data\input\galapagos_field_data\GridsDataFrame.csv')
        self.TakeOffPercentages = None
        self.TakeOffAbsolute = None
        
        self.BeachingPercentages = None
        self.BeachingAbsolute = None
        
    # this method computes the take-off locations for each coast grid. 
    def take_off_locations(self):
        
        #To compute this, we need the grids dataframe and the transition matrix
        gridsdataframe = self.GridDataFrame
        transitionmatrix = self.TransitionMatrix
        
        #an array which will be filled with the absolute number of particles that take-off from this island an beaches on another island. 
        abs_numb_to_other_island = np.zeros((len(gridsdataframe)))

        #loop over all release grids
        for release_grid in range(len(transitionmatrix[:,1])):
            #check for each release grid to which island it belongs
            release_island = gridsdataframe['island_number'][release_grid]
            
            #loop over all beach grids
            for beach_grid in range(len(transitionmatrix[1,:])):
                
                #check for each beach grid to which island it belongs
                beach_island = gridsdataframe['island_number'][beach_grid]
                print('beach island is:', beach_island)
        
                #if the release island is NOT the same as the beach island we take these particles into account
                if release_island != beach_island:
#                    print('test')
                    #check how many particles have been trasported from the release grid to the beach grid
                    number_particles = transitionmatrix[release_grid, beach_grid]
                    #Add this number of partciles to the entry of abs_numb_to_other_island array that corresponds to this release grid
                    abs_numb_to_other_island[release_grid] += number_particles
                
        #compute the number of particles that are released from the release grid have beached somewhere. So doesn't matter if they beached on
        #the same island they were released from or on another islands
        beached_part = np.sum(transitionmatrix, axis = 1)
        
        #Here we compute the take-off percentage, which represent percentage of particles that have beached on another island relative to
        #the total number of particles that have beached.
        self.TakeOffPercentages = (abs_numb_to_other_island/beached_part)*100
        
        self.TakeOffAbsolute = abs_numb_to_other_island
        
        
    #here compute the beach locations. Everything we do is exaclty the same except, that we transpose the 
    #transition matrix
    def beach_locations(self):
        
        gridsdataframe = self.GridDataFrame
        transitionmatrix = self.TransitionMatrix
        
        transitionmatrix_transpose = transitionmatrix.transpose()
        
        abs_numb_to_other_island = np.zeros((len(gridsdataframe)))

        for release_grid in range(len(transitionmatrix_transpose[:,1])):
            release_island = gridsdataframe['island_number'][release_grid]
        #    print('release island is:', release_island)
            for beach_grid in range(len(transitionmatrix_transpose[1,:])):
                
                beach_island = gridsdataframe['island_number'][beach_grid]
                print('beach island is:', beach_island)
        
                if release_island != beach_island:
                    print('test')
                    
                    number_particles = transitionmatrix_transpose[release_grid, beach_grid]
                    abs_numb_to_other_island[release_grid] += number_particles
                
        beached_part = np.sum(transitionmatrix_transpose, axis = 1)
        
        self.BeachingPercentages = (abs_numb_to_other_island/beached_part)*100
        self.BeachingAbsolute = abs_numb_to_other_island
        
    
    #here we plot the beaching and take-off locations. 
    def plot(self, savename_fig, beach = None, takeoff = None):
        
        if takeoff == True:
        
            scaling = 0.3
            NetworkDataFrame = pd.DataFrame({'GridNumber':self.GridDataFrame['Grid_Number'], 'longitude':self.GridDataFrame['min_lon'], 'latitude':self.GridDataFrame['min_lat'],
                                             'abs_number_takeoff': self.TakeOffAbsolute , 'perc_number_takeoff':self.TakeOffPercentages})
    
            plt.figure(figsize = (10,10))
            ax = plt.axes(projection = ccrs.PlateCarree())
            ax.set_extent([-92, -89, -1.5, 0.75], ccrs.PlateCarree())         
            ax.add_feature(cfeature.LAND)
            ax.add_feature(cfeature.OCEAN)
            ax.add_feature(cfeature.COASTLINE)
            grid_lines = ax.gridlines(draw_labels=True, linestyle = '-.', color = 'black')
            grid_lines.xformatter = LONGITUDE_FORMATTER
            grid_lines.yformatter = LATITUDE_FORMATTER
            sc = plt.scatter(NetworkDataFrame['longitude'], NetworkDataFrame['latitude'], c = NetworkDataFrame['perc_number_takeoff'], s = NetworkDataFrame['abs_number_takeoff']*scaling)
            plt.set_cmap('Reds')
            plt.colorbar(orientation="horizontal").set_label('Percentage particles beached on another island (%)', fontsize = 15)
            kw = dict(prop="sizes", num=5, color=sc.cmap(0.7))
            plt.legend(*sc.legend_elements(**kw), title = 'Absolute number of particles', loc = 'upper right', fontsize = 12)
            plt.savefig(parent_dir + '\Figures\Take_off_locations\Take_off_' + savename_fig)
            
            
        elif beach == True:
            
            scaling = 0.3
            NetworkDataFrame = pd.DataFrame({'GridNumber':self.GridDataFrame['Grid_Number'], 'longitude':self.GridDataFrame['min_lon'], 'latitude':self.GridDataFrame['min_lat'],
                                             'abs_number_beach': self.BeachingAbsolute , 'perc_number_beach':self.BeachingPercentages})
    
            plt.figure(figsize = (10,10))
            ax = plt.axes(projection = ccrs.PlateCarree())
            ax.set_extent([-92, -89, -1.5, 0.75], ccrs.PlateCarree())         
            ax.add_feature(cfeature.LAND)
            ax.add_feature(cfeature.OCEAN)
            ax.add_feature(cfeature.COASTLINE)
            grid_lines = ax.gridlines(draw_labels=True, linestyle = '-.', color = 'black')
            grid_lines.xformatter = LONGITUDE_FORMATTER
            grid_lines.yformatter = LATITUDE_FORMATTER
            sc = plt.scatter(NetworkDataFrame['longitude'], NetworkDataFrame['latitude'], c = NetworkDataFrame['perc_number_beach'], s = NetworkDataFrame['abs_number_beach']*scaling)
            plt.set_cmap('Reds')
            plt.colorbar(orientation="horizontal").set_label('Percentage particles beached on another island (%)', fontsize = 15)
            kw = dict(prop="sizes", num=5, color=sc.cmap(0.7))
            plt.legend(*sc.legend_elements(**kw), title = 'Absolute number of particles', loc = 'upper right', fontsize = 12)
            plt.savefig(parent_dir + '\Figures\Beaching_locations\Beaching_' + savename_fig)
            
            
                
#
#        
#
#        
        
    