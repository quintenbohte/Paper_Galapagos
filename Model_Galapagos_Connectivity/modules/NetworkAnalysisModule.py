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
from pathlib import Path
import os
parent_dir = str(Path(os.path.abspath(__file__)).parents[1])


#This class contains the methods which executes the network analysis.
class NetworkAnalysis:
    
    def __init__(self, TransitionMatrix, GridsDataFrame):

        self.TransitionMatrix = TransitionMatrix
        self.Graph = None
        self.GridsDataFrame = GridsDataFrame
        self.NetworkDataFrame = None
        self.BetweennessCentralities = None
        self.BetweennessCentralitiesNorm = None
        
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
            #loop over all nodes in the connstructed graph, these will represent the release nodes
            for release_node in ConstructedGraph.nodes():

                #loop over all nodes in the constructed graph, these will represent the beach nodes
                for beach_node in ConstructedGraph.nodes():
                    
                    #We loop over all release nodes and all beach nodes. Now we want to know if there is an edge between them. For this, the
                    #transition matrix is used. Using the number of the release grid and the number of the beach grid we can find the corresponding entry in the
                    #transition matrix, the value of this entry is the weight of the edge between the release and beach node.
                    weight = TransitionMatrix[release_node, beach_node]
                    #compute the weight log, which is needed for computing the most likely path
                    weightlog = -np.log(TransitionMatrix[release_node, beach_node])
                    
                    #If the value of the transition matrix is <0, no particles are transported from this release node to this
                    #beach node, so no edge will be added. However, if weight > 0:
                    if weight > 0:
                        #we add an adge, and the weight of this edge and the log weight of this edge. 
                        ConstructedGraph.add_edge(release_node,
                                       beach_node,
                                       weightpath = weight,
                                       weightlogpath = weightlog)
        #here we do the exact same as above, but now with all edges reversed. Note that to do this we only have to transpose our transition matrix.
        if ReverseDirections == True:
            for release_node in ConstructedGraph.nodes():

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


    #This function will compute the betweenness centrality of each node using the most likely path    
    def BetwneennesCentrality(self):
        #the graph we will use
        ConstructedGraph = self.Graph
        
        #the total number of nodes
        number_of_nodes = ConstructedGraph.number_of_nodes()
        
        #the list which will be filled with all the most likely paths between all nodes
        PathList = []
        #the list which will be filled with all the beach/target nodes
        targetList = []
        #the list which will be filled with all the release/source nodes
        sourceList = []
        
        #loop over all nodes (representing the target nodes)
        for target in range(number_of_nodes):
            
            #loop over all nodes (representing the source nodes)
            for source in range(number_of_nodes):
                
                #we use the try and except command, since if there is no path between two nodes it will give an error. With the try command
                #it will then just go to the next one
                try:
                    #compute the ShortestPath with the dijkstra algorithm.
                    ShortestPathDijkstraLog = nx.dijkstra_path(ConstructedGraph, target,source, weight = 'weightlogpath')
                    #fill all lists we defined above
                    targetList.append(target)
                    sourceList.append(source)
                    #Note, this list will contain all most likely paths between all different nodes
                    PathList.append(ShortestPathDijkstraLog)
                    
                except:
                    pass
            
        #this will be the list with the betweenness centrality for each node
        BetwCentrList = []
        #loop over all nodes
        for grid in range(number_of_nodes):
            #this is the count for the number of shortest path that cross a node
            count = 0
            #we loop over all shortest paths, that we found in our graph
            for path in PathList:
                
                #Now we check for every grid/node in how many shortest paths it appears
                if grid in path:
                    #for every time the grid/node appears in a shorest path, we add 1 to the count
                    count += 1
                    
            #The betweenness centrality is the percentage of shorest paths that cross a node. So we divide the count by the
            #total number of shortest paths. 
            BetwCentr = count/len(PathList)
            
            BetwCentrList.append(BetwCentr)
            
            
        self.BetweennessCentralities = BetwCentrList
        data = self.BetweennessCentralities
        self.BetweennessCentralitiesNorm =  (data - np.min(data)) / (np.max(data) - np.min(data))
#        return BetwCentrList

    #this method computes the normal betweenness centrality, using the networkx module.
    def BetwneennesCentralityNetworkx(self):
         ConstructedGraph = self.Graph
         self.BetweennessCentralities = nx.betweenness_centrality(ConstructedGraph)
         data = self.BetweennessCentralities
         self.BetweennessCentralitiesNorm =  (data - np.min(data)) / (np.max(data) - np.min(data))

    

    #This method make a a dataframe, where each row represents a node. The column values contain the lon, lat, grid number, and betweenness centrality of the node
    #this dataframe is very convenient, since it contain all usefull information of the network analysis and it is usefull for plotting. 
    def ConstructNetworkDataFrame(self):

        BetwCentrList = self.BetweennessCentralities
        BetwCentrListNorm = self.BetweennessCentralitiesNorm
        GridsDataFrame = self.GridsDataFrame

        NetworkDataFrame = pd.DataFrame({'GridNumber':GridsDataFrame['Grid_Number'], 'longitude':GridsDataFrame['min_lon'], 'latitude':GridsDataFrame['min_lat'],
                                 'Betweennes_Centrality': BetwCentrList, 'Betweennes_Centrality_norm':BetwCentrListNorm})
    
        self.NetworkDataFrame = NetworkDataFrame


    #This method plots the betweenness centrality on a map of the galapagos. 
    def PlotCentrality(self, scaling, savename, title, log = False):
        
        NetworkDataFrame = self.NetworkDataFrame
        
        if log == False:
            centrality = NetworkDataFrame['Betweennes_Centrality']
            plt.figure(figsize = (10,10))
            ax = plt.axes(projection = ccrs.PlateCarree())
            ax.set_extent([-92, -89, -1.5, 0.75], ccrs.PlateCarree())         
            ax.add_feature(cfeature.LAND)
            ax.add_feature(cfeature.OCEAN)
            ax.add_feature(cfeature.COASTLINE)
            grid_lines = ax.gridlines(draw_labels=True, linestyle = '-.', color = 'black')
            grid_lines.xformatter = LONGITUDE_FORMATTER
            grid_lines.yformatter = LATITUDE_FORMATTER
            plt.scatter(NetworkDataFrame['longitude'], NetworkDataFrame['latitude'], c = centrality, s = centrality*scaling)
    #        plt.scatter(NetworkDataFrame['longitude'], NetworkDataFrame['latitude'], c = centrality, s = centrality)
            plt.set_cmap('Reds')
            plt.title(str(title), fontsize = 20)
            plt.colorbar(orientation="horizontal").set_label('Betweenness Centrality (percentage most likely paths crossed)', fontsize = 15)
            plt.savefig(parent_dir + '\Figures\Betweenness_Centrality\Betw_cen_' + savename + '_' + title)
     
            
            
        elif log == True:
            centrality = NetworkDataFrame['Betweennes_Centrality_norm']
        
            plt.figure(figsize = (10,10))
            ax = plt.axes(projection = ccrs.PlateCarree())
            ax.set_extent([-92, -89, -1.5, 0.75], ccrs.PlateCarree())         
            ax.add_feature(cfeature.LAND)
            ax.add_feature(cfeature.OCEAN)
            ax.add_feature(cfeature.COASTLINE)
            grid_lines = ax.gridlines(draw_labels=True, linestyle = '-.', color = 'black')
            grid_lines.xformatter = LONGITUDE_FORMATTER
            grid_lines.yformatter = LATITUDE_FORMATTER
            plt.scatter(NetworkDataFrame['longitude'], NetworkDataFrame['latitude'], c = -1*np.log(centrality), s = scaling*(-1*np.log(centrality)))
            plt.set_cmap('Reds')
    #        plt.scatter(NetworkDataFrame['longitude'], NetworkDataFrame['latitude'], c = centrality, s = centrality)
            plt.title(str(title), fontsize = 20)
            plt.colorbar(orientation="horizontal").set_label('Normalized Betweenness Centrality', fontsize = 15)
            plt.savefig(parent_dir + '\Figures\Betweenness_Centrality\Betw_cen_' + savename + '_' + title)
            
            
        return
    
    

    




























