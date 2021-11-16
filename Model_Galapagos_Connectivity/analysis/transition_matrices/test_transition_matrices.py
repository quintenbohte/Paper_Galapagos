# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 10:32:13 2021

@author: quint
"""

import os
from pathlib import Path
import sys
parent_dir = str(Path(os.path.abspath(__file__)).parents[2])
sys.path.append(parent_dir + '\modules')
#sys.path.append(r"C:\Users\quint\Documents\Quinten_studie\Publicatie\Modules")

import numpy as np
import pandas as pd
from SimulationModule import Simulation as SimM
#from TrajectoryPlotsModule import TrajectoryPlots as TP
from SimpleFunctions import Functions as func
from TransitionMatrixModuleLoc import ConstructTransitionMatrixOnlyTrajData as TM
from NetworkAnalysisModule import NetworkAnalysis as NW
from TrajectoryPlotsModuleLoc import TrajectoryPlots as TP
import matplotlib.pyplot as plt

savename = 'transition_matrix_test'
parent_dir = str(Path(os.path.abspath(__file__)).parents[2])
path_to_dict = parent_dir + '\data\output\simulations\Simulation_' + savename + '.dictionary'
GridsDataFrame = pd.read_csv(parent_dir + '\data\input\galapagos_field_data\GridsDataFrame.csv')
traj_data_dict = func.load_dict(path_to_dict)
domain = [-92, -88, -2, 2]

Plot_obj = TP(savename = savename , domain = domain)
Plot_obj.plot_trajectories()









#%%
#
#delta_time = 1
#beachingtimescale = 5
#
#tm_obj = TM(traj_data_dict, delta_time, savename, beachingtimescale)
#tm_obj.construct_matrix()
#transition_matrix = tm_obj.transition_matrix
#
#




