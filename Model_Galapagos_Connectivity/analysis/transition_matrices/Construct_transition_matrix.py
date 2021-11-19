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


savename = '2008_2012'
path_to_dict = '/data/oceanparcels/output_data/data_Quinten/Paper_Galapagos/Data/Simulations/Simulation2008-2012.dictionary'

traj_data_dict = func.load_dict(path_to_dict)

Tm_obj_wob = TraM(traj_data_dict, delta_time = 1, savename = savename, beaching_timescale = 5)
Tm_obj_wob.construct_matrix()
tm = Tm_obj_wob.transition_matrix

np.save('/data/oceanparcels/output_data/data_Quinten/Paper_Galapagos/Data/Simulations/TransitionMatrix2008_2012', tm)



