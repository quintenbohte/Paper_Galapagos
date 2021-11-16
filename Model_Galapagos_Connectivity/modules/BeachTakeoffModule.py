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



class BeachTakeOffAnalysis():
    
    
    def __init__(self, TransitionMatrix):

        self.TransitionMatrix = TransitionMatrix
        self.GridDataFrame = pd.read_csv(parent_dir + '\data\input\galapagos_field_data\GridsDataFrame.csv')
        self.TakeOffPercentages = None
        self.TakeOffAbsolute = None
        
        self.BeachingPercentages = None
        self.BeachingAbsolute = None
        
    def take_off_locations(self):
        
        gridsdataframe = self.GridDataFrame
        transitionmatrix = self.TransitionMatrix
        
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
                
        beached_part = np.sum(transitionmatrix, axis = 1)
        
        self.TakeOffPercentages = (abs_numb_to_other_island/beached_part)*100
        self.TakeOffAbsolute = abs_numb_to_other_island
        
        
        
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
        
    