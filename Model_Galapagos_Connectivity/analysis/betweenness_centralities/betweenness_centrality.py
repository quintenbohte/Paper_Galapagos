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

savename = 'long_run'
transitionmatrix = np.load(parent_dir + '\data\output\Transition_matrices\Tm_' +savename + '.npy')
gridsdataframe = pd.read_csv(parent_dir + '\data\input\galapagos_field_data\GridsDataFrame.csv')
np.fill_diagonal(transitionmatrix,0)


scaling = 12000
savename_fig = 'betw_centr_long_run'
NetworkAnalysisObject = Nw(transitionmatrix, gridsdataframe)
NetworkAnalysisObject.CreateDirectedGraph()
NetworkAnalysisObject.BetwneennesCentrality()
NetworkAnalysisObject.ConstructNetworkDataFrame()
NetworkAnalysisObject.PlotCentrality(scaling = scaling, savename=savename_fig, title = 'betweenness centrality', log = False)


















