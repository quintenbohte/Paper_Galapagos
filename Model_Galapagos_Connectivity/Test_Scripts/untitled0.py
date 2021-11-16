# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 10:43:53 2021

@author: quint
"""


from parcels import Field, FieldSet, JITParticle, ScipyParticle 
from parcels import ParticleFile, ParticleSet, Variable, VectorField, ErrorCode
from parcels.tools.converters import GeographicPolar 
from datetime import timedelta as delta
from os import path
from glob import glob
import numpy as np
import dask
import os

import math
import xarray as xr
#from parcels import AdvectionRK4
from netCDF4 import Dataset
import warnings
import matplotlib.pyplot as plt
import pickle
warnings.simplefilter('ignore', category=xr.SerializationWarning)
from operator import attrgetter
from pathlib import Path
import sys


parent_dir = str(Path(os.path.abspath(__file__)).parents[1])


sys.path.append(r"C:\Users\quint\Documents\Quinten_studie\Publicatie\Modules")



savename = 'test_interpolation_delete_part_3'


Sim_data = xr.open_dataset(parent_dir + '\data\output\simulations\Simulation_' + savename + '.nc')
lon_traj = Sim_data['lon'].data

#%%

#savename = 'test_interpolation_delete_part'
#
#Sim_data_rem = xr.open_dataset(parent_dir + '\data\output\simulations\Simulation_' + savename + '.nc')
#lon_traj_rem = Sim_data_rem['lon'].data


def replace_removed_part(sim_data, removed_part_number):
    rows = len(sim_data[:,0])
    columns = len(sim_data[1,:])
    
    sim_data_new = np.zeros((rows+1,columns))
    
    sim_data_new[0:removed_part_number-1,:] = sim_data[0:removed_part_number-1,:]
    sim_data_new[removed_part_number+1:,:] = sim_data[removed_part_number:,:]
    
    return sim_data_new
    

sim_data_new = replace_removed_part(lon_traj, 23000)
    
    
    
subset = lon_traj[22500:23500,:]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




