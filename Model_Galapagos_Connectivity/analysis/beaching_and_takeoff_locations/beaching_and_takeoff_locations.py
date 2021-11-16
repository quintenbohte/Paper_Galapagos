# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 15:01:51 2021

@author: quint
"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import os

parent_dir = str(Path(os.path.abspath(__file__)).parents[2])
sys.path.append(parent_dir + '\modules')
from BeachTakeoffModule import BeachTakeOffAnalysis as BTA

#########################################DEFINE RELATIVE PATH#######################################

parent_dir = str(Path(os.path.abspath(__file__)).parents[2])

########################################## LOAD DATA ######################################

savename = 'long_run'
transitionmatrix = np.load(parent_dir + '\data\output\Transition_matrices\Tm_' +savename + '.npy')
gridsdataframe = pd.read_csv(parent_dir + '\data\input\galapagos_field_data\GridsDataFrame.csv')

########################################## DETERMINE BEACHING AND TAKEOFF ###########################################

BeachTakeoffObject = BTA(transitionmatrix)
BeachTakeoffObject.take_off_locations()
BeachTakeoffObject.beach_locations()

BeachTakeoffObject.plot(savename_fig = 'long_run', takeoff = True)
BeachTakeoffObject.plot(savename_fig = 'long_run', beach = True)

















