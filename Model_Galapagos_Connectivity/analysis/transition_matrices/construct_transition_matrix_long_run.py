# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 16:19:34 2021

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

def Convert_to_single_particles(data):

    data = data
    output_frequency = 1     #hours
    repeated_days = 60
    advection_duration = int(repeated_days*(24/output_frequency))+1 
    release_locations = 393 
#    length_simulation = 1450 #days
    advected_timesteps = len(data[1,:])
    number_particles = math.ceil(advected_timesteps/advection_duration)
    
    
    output_trajectory = np.zeros(((len(data[:,1])*number_particles) ,advection_duration))
    
    for trajectory_number in range(len(data)):
        
        trajectory = data[trajectory_number,:]
    
        remaining_time_steps = len(trajectory)%advection_duration
    
        full_particle_trajectories = trajectory[0:len(trajectory)-remaining_time_steps]
        
        seperate_particles = np.split(full_particle_trajectories,number_particles-1)
        last_particle = trajectory[len(full_particle_trajectories):(len(full_particle_trajectories)+remaining_time_steps)]
    
        particle_set = int(math.floor(trajectory_number/release_locations))
    
        for particle_number in range((number_particles)): 
            
            if particle_number < number_particles-1:
            
                print('part_set:',particle_set)
                indices_current_particle_set = particle_set * (release_locations*number_particles)
                print('current_ind:',indices_current_particle_set)
                index = indices_current_particle_set + ((trajectory_number%release_locations) + (release_locations*particle_number))
                print('index:',index)
                
                output_trajectory[index,:] = seperate_particles[particle_number]
            
            else:
                print('part_set:',particle_set)
                indices_current_particle_set = particle_set * (release_locations*number_particles)
                print('current_ind:',indices_current_particle_set)
                index = indices_current_particle_set + ((trajectory_number%release_locations) + (release_locations*particle_number))
                print('index:',index)
                
                output_trajectory[index,0:217] = last_particle
            
    
    return output_trajectory


Sim_data = xr.open_dataset(parent_dir + '\data\output\simulations\Simulation2008_2010_2011_2012.nc')

lon = Sim_data['lon'].data
lat = Sim_data['lat'].data
coastgrid = Sim_data['coastcell'].data
beached = Sim_data['beached'].data
grid_number = Sim_data['gridnumber'].data

coastgrids_data = Convert_to_single_particles(coastgrid)
grid_number_data = Convert_to_single_particles(grid_number)
traj_data_dict = {'coastcells':coastgrids_data, 'gridnumber':grid_number_data}

savename =  'long_run'
Tm_obj_wob = TraM(traj_data_dict, delta_time = 1, savename = savename, beaching_timescale = 5)
Tm_obj_wob.construct_matrix()
tm = Tm_obj_wob.transition_matrix












