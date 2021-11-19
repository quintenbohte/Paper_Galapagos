# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 12:25:06 2021

@author: quint
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
import os
parent_dir = str(Path(os.path.abspath(__file__)).parents[2])
sys.path.append(parent_dir + '\modules')
from TransitionMatrixModuleLoc import ConstructTransitionMatrixOnlyTrajData as TraM
import math
import xarray as xr
from SimpleFunctions import Functions as func


#define the path to the dictionary with the simulation data.
#the savename should be the same savename as you used in for the simulation
savename = 'test'
path_to_dict = parent_dir + '\data\output\simulations\Simulation_' + savename +  '.dictionary'

#load the dictionary
traj_data_dict = func.load_dict(path_to_dict)

#the savename is the name the of the file after we save the transitionmatrix

#make the transiton matrix object, Tm_obj
Tm_obj= TraM(traj_data_dict, delta_time = 1, savename = savename, beaching_timescale = 5)

#use the construct_matrix method to make the transition matrix. 
#See TransitionMatrixModuleLoc, in folder modules, for the full code of this method.
#Note, we use the ConstructTransitionMatrixOnlyTrajData class. NOT the ConstructTransitionMatrix class.

Tm_obj.construct_matrix()

#get the transition matrix from the Tm_object
tm = Tm_obj.transition_matrix

#save the transition matrix
np.save(parent_dir + '\data\output\Transition_matrices\Tm_' + savename, tm)




