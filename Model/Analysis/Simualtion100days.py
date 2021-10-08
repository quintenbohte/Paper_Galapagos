# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 09:46:46 2021

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
import matplotlib.pyplot as plt

#%%

length_simulation = 2 #unit: days (for how long do we deploy particles)
advection_duration = 60 #unit: days (how long does one particle advect in the fields)
output_frequency = 1     #unit: hours
repeatdt = 24        #unit: hours
deltatime = 1           #dt in hours
domain = [-92, -88, -2, 2]
beaching_timescale = 3

data_in = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input'
path_to_velocity = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input\RGEMS_2010_Surf.nc'
path_to_grid = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input\RGEMS3_Surf_grid.nc'
savename = '1year'
###############################################################
#                    RUN SIMULATION                           #
###############################################################


sim1year = SimM(length_simulation = length_simulation, advection_duration = advection_duration,
                              output_frequency = output_frequency, repeatdt = repeatdt, 
                              deltatime = deltatime, savename = savename,
                              domain = domain,
                              path_to_velocity = path_to_velocity,
                              path_to_grid = path_to_grid, data_in=data_in)

sim1year.run_simulation()


#%%


