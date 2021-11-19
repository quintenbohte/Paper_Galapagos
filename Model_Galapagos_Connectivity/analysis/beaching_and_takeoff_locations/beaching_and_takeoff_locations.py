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

#savename is the name we gave as input when the transition matrix was created
savename = 'test'
#load the transition matrix
transitionmatrix = np.load(parent_dir + '\data\output\Transition_matrices\Tm_' +savename + '.npy')
#load the gridsdataframe
gridsdataframe = pd.read_csv(parent_dir + '\data\input\galapagos_field_data\GridsDataFrame.csv')

########################################## DETERMINE BEACHING AND TAKEOFF ###########################################

#make beachtakeoffobject using the BeachTakeOff Analysis class (BTA). Only takes the transition matrix as input
BeachTakeoffObject = BTA(transitionmatrix)

#compute the take off locations, using the take_off_locations method in the BTA
BeachTakeoffObject.take_off_locations()

#compute the beach_locations, using the the beach_locations method
BeachTakeoffObject.beach_locations()

#plot the takeoff and beach locations. 
savename_fig = 'test'
BeachTakeoffObject.plot(savename_fig = savename_fig, takeoff = True)
BeachTakeoffObject.plot(savename_fig = savename_fig, beach = True)

















