# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 16:25:26 2021

@author: quint
"""

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

class NetworkAnalysis:
    
    def __init__(self, TransitionMatrix, GridsDataFrame):

        self.TransitionMatrix = TransitionMatrix
        self.Graph = None
        self.GridsDataFrame = GridsDataFrame
        self.NetworkDataFrame = None
        self.BetweennessCentralities = None
        
    'ClassMethods'
    
    def CreateDirectedGraph(self, ReverseDirections = False, DrawNetwork = False):
        
        GridsDataFrame = self.GridsDataFrame
        TransitionMatrix = self.TransitionMatrix
        '''
        This function will create a weighted directed graph from a transition matrix
        '''
        '''
        load the dataframe with the location (loongitude and latitude) of the nodes
        '''
        
        '''
        This loop creates two lists with the locations and nodenames (numbers in our case).
        '''
        nodenames = []
        locations = []
        for i in range(len(GridsDataFrame)):
            lon =  GridsDataFrame['min_lon'].iloc[i]
            lat =  GridsDataFrame['min_lat'].iloc[i]
            location = (lon,lat)
            locations.append(location)
            nodenames.append(i)
        '''
        create an empty graph
        '''
        ConstructedGraph = nx.DiGraph()
        
        '''
        loop over all grids and make a node for every coast grid. We can add attributes to the nodes. In this case we added the attribute location. (and attribute shit).
        you can take a look at the attributes by typing Gd.nodes('name of attribute'). Example: Gd.nodes('locations')
        '''
        for i in range(len(GridsDataFrame)):
            nodename = i
            loc = locations[i]
            ConstructedGraph.add_node(nodename, pos=loc)
        '''
        here we substract the location attribute for plotting the graph
        '''
#        pos = nx.get_node_attributes(ConstructedGraph, 'pos')
        '''
        add the edges to the graph. We add to attributes to the edge:
            - normal weight
            - and the log of the normal weight, this one is used for calculating the most likely path
        '''
        
        if ReverseDirections == False:
            for release_node in ConstructedGraph.nodes():
                print(release_node)
                for beach_node in ConstructedGraph.nodes():
                    weight = TransitionMatrix[release_node, beach_node]
                    weightlog = -np.log(TransitionMatrix[release_node, beach_node])
                    if weight > 0:
                        ConstructedGraph.add_edge(release_node,
                                       beach_node,
                                       weightpath = weight,
                                       weightlogpath = weightlog)
                        
        if ReverseDirections == True:
            for release_node in ConstructedGraph.nodes():
                print(release_node)
                for beach_node in ConstructedGraph.nodes():
                    weight = TransitionMatrix.T[release_node, beach_node]
                    weightlog = -np.log(TransitionMatrix.T[release_node, beach_node])
                    if weight > 0:
                        ConstructedGraph.add_edge(release_node,
                                       beach_node,
                                       weightpath = weight,
                                       weightlogpath = weightlog)
            
    
        self.Graph = ConstructedGraph
#        return ConstructedGraph
        
        
    
    
    def DrawGraph(self, ConstructedGraph):
        
        
        GridsDataFrame = self.GridsDataFrame

        
        nodenames = []
        locations = []
        
        for i in range(len(GridsDataFrame)):
            lon =  GridsDataFrame['min_lon'].iloc[i]
            lat =  GridsDataFrame['min_lat'].iloc[i]
            location = (lon,lat)
            locations.append(location)
            nodenames.append(i)
            
                
        for i in range(len(GridsDataFrame)):
            nodename = i
            loc = locations[i]
            ConstructedGraph.add_node(nodename, pos=loc)
        
        pos = nx.get_node_attributes(ConstructedGraph, 'pos')
        
        
        plt.figure(figsize = (10,10))
        ax = plt.axes(projection = ccrs.PlateCarree()) 
        ax.set_extent([-92, -89, -1.5, 0.5], ccrs.PlateCarree())         
        ax.add_feature(cfeature.LAND, facecolor = 'green')
        ax.add_feature(cfeature.OCEAN)
        ax.add_feature(cfeature.COASTLINE)
        ax.set_title('Directed Network', fontsize =20)    
        nx.draw_networkx_edges(ConstructedGraph, pos)
    
        return 
    
    def BetwneennesCentrality(self):

        ConstructedGraph = self.Graph
        
        number_of_nodes = ConstructedGraph.number_of_nodes()
        
        PathList = []
        targetList = []
        sourceList = []
        for target in range(number_of_nodes):
            print(target)
            for source in range(number_of_nodes):
                
                try:
                    ShortestPathDijkstraLog = nx.dijkstra_path(ConstructedGraph, target,source, weight = 'weightlogpath')
                    targetList.append(target)
                    sourceList.append(source)
                    PathList.append(ShortestPathDijkstraLog)
                    
                except:
                    pass
            
            
        BetwCentrList = []
        for grid in range(number_of_nodes):
            print(grid)
            count = 0
            for path in PathList:
        
                if grid in path:
                    count += 1
                    
            BetwCentr = count/len(PathList)
            
            BetwCentrList.append(BetwCentr)
            
            
        self.BetweennessCentralities = BetwCentrList
        
#        return BetwCentrList


    def ConstructNetworkDataFrame(self):

        BetwCentrList = self.BetweennessCentralities
        GridsDataFrame = self.GridsDataFrame

        NetworkDataFrame = pd.DataFrame({'GridNumber':GridsDataFrame['Grid_Number'], 'longitude':GridsDataFrame['min_lon'], 'latitude':GridsDataFrame['min_lat'],
                                 'Betweennes_Centrality': BetwCentrList})
        self.NetworkDataFrame = NetworkDataFrame



    def PlotCentrality(self, scaling, title):
        
        
        NetworkDataFrame = self.NetworkDataFrame
        
        plt.figure(figsize = (10,10))
        ax = plt.axes(projection = ccrs.PlateCarree())
        ax.set_extent([-92, -89, -1.5, 0.75], ccrs.PlateCarree())         
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.OCEAN)
        ax.add_feature(cfeature.COASTLINE)
        grid_lines = ax.gridlines(draw_labels=True, linestyle = '-.', color = 'black')
        grid_lines.xformatter = LONGITUDE_FORMATTER
        grid_lines.yformatter = LATITUDE_FORMATTER
        plt.scatter(NetworkDataFrame['longitude'], NetworkDataFrame['latitude'], c = NetworkDataFrame['Betweennes_Centrality'], s = NetworkDataFrame['Betweennes_Centrality']*scaling)
        plt.title(str(title))
        
        return
            
    































