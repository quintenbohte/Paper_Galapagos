# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 12:31:00 2021

@author: quint
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
import os
parent_dir = str(Path(os.path.abspath(__file__)).parents[2])
sys.path.append(parent_dir + '\modules')

from NetworkAnalysisModule import NetworkAnalysis as Nw

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt


#########################################DEFINE RELATIVE PATH#######################################

parent_dir = str(Path(os.path.abspath(__file__)).parents[2])

########################################## LOAD DATA ######################################

#savename is the name we gave to the transition matrix when we created it. 
savename = 'long_run'
transitionmatrix = np.load(parent_dir + '\data\output\Transition_matrices\Tm_' +savename + '.npy')

#load gridsdataframe, this includes all information about the coastgrids.
gridsdataframe = pd.read_csv(parent_dir + '\data\input\galapagos_field_data\GridsDataFrame.csv')

#make the diagonal zero, so we do not include the particles which beaches at their grid where they were released from
np.fill_diagonal(transitionmatrix,0)

#scaling only makes the scatters larger so they are more clearlt visible in the plot. 
scaling = 12000
savename_fig = 'betw_centr_long_run'

#create NetworkAnalysisObject using the NetworkAnalysis class. It takes as input the transition matrix and the gridsdataframe. This object includes all methods to execute our
#network analysis
NetworkAnalysisObject = Nw(transitionmatrix, gridsdataframe)

#first we create a directed graph from the transition matrix
NetworkAnalysisObject.CreateDirectedGraph()

#compute the betweenness centrality, using the directed graph
NetworkAnalysisObject.BetwneennesCentrality()

#construct networkdataframe, this is needed to make the plots
NetworkAnalysisObject.ConstructNetworkDataFrame()

#plot the betweenness centrality on a map
NetworkAnalysisObject.PlotCentrality(scaling = scaling, savename=savename_fig, title = 'betweenness centrality', log = False)


















