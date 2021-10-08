# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 11:28:58 2021

@author: quint
"""

import sys
sys.path.append(r"C:\Users\quint\Documents\Quinten_studie\Publicatie\Modules")
import numpy as np
import pandas as pd
from SimulationModule import Simulation as SimM
from TrajectoryPlotsModule import TrajectoryPlots as TP
from SimpleFunctions import Functions as func
from TransitionMatrixModule import ConstructTransitionMatrix as TM
from NetworkAnalysisModule import NetworkAnalysis as NW
import matplotlib.pyplot as plt


path_to_dict = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Simulations\dict100Days_sim.dictionary'
GridsDataFrame = pd.read_csv(r'C:\Users\quint\Documents\Quinten_studie\Master_Thesis\Data\Output_Data\ReleaseLocations\GridsDataFrame.csv')
sim_data = func.load_dict(path_to_dict)
beachingtimescale = 5

tm_obj = TM(sim_data,beachingtimescale)
tm_obj.construct_matrix()
transition_matrix = tm_obj.transition_matrix

#%%
nwobj = NW(transition_matrix,GridsDataFrame)
nwobj.CreateDirectedGraph()
nwobj.BetwneennesCentrality()
betw_centr = nwobj.BetweennessCentralities



