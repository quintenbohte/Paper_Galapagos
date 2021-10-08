# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 14:19:01 2021

@author: quint
"""

import sys
sys.path.append(r"C:\Users\quint\Documents\Quinten_studie\Publicatie\Modules")
import numpy as np
import pandas as pd
from SimulationModule import Simulation as SimM
from TrajectoryPlotsModule import TrajectoryPlots as TP
from TransitionMatrixModule import ConstructTransitionMatrix as TraM
import matplotlib.pyplot as plt

###############################################################
#                    SIMULATION VARIABLES                     #
###############################################################


length_simulation = 20 #unit: days (for how long do we deploy particles)
advection_duration = 60 #unit: days (how long does one particle advect in the fields)
output_frequency = 1     #unit: hours
repeatdt = 24        #unit: hours
deltatime = 1           #dt in hours
domain = [-92, -88, -2, 2]
beaching_timescale = 3

data_in = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input'
path_to_velocity = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input\RGEMS_2010_Surf.nc'
path_to_grid = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input\RGEMS3_Surf_grid.nc'
#%%
###############################################################
#                    RUN SIMULATION                           #
###############################################################

bordercurrentlist = [0, 0.0000003, 0.0000006,0.0000009,0.0000012,0.0000015,0.0000017,0.0000019,0.0000021,0.0000023]


for bordercurrent in bordercurrentlist:
    savename = str(bordercurrent)
    sim_wbc = SimM(length_simulation = length_simulation, advection_duration = advection_duration,
                                  output_frequency = output_frequency, repeatdt = repeatdt, 
                                  deltatime = deltatime, savename = savename,
                                  domain = domain,
                                  path_to_velocity = path_to_velocity,
                                  path_to_grid = path_to_grid, data_in=data_in, bordercurrent = bordercurrent)
    sim_wbc.run_simulation()
    beached_particles_per_border_current.append(sim_wbc.beached_particles)
    TP(sim_wbc).plot_trajectories()
    

plt.figure()
plt.plot(bordercurrentlist, beached_particles_per_border_current)
plt.axvline(x=0.0000012, color = 'red', label='line at x = {}'.format(0.0000012))
plt.ylabel('percentage particles on land')
plt.xlabel('value border current')
plt.legend()

#%%
################################################


bordercurrent = 0.00000000
savename = 'wob'
sim_wobc = SimM(length_simulation = length_simulation, advection_duration = advection_duration,
                              output_frequency = output_frequency, repeatdt = repeatdt, 
                              deltatime = deltatime, savename = savename,
                              domain = domain,
                              path_to_velocity = path_to_velocity,
                              path_to_grid = path_to_grid, data_in=data_in, bordercurrent = bordercurrent)
sim_wobc.run_simulation()
Tm_obj_wob = TraM(sim_wobc, beaching_timescale)


#%%
length = [320]
#beached_particles_list = []
bordercurrent = 0.00000000
for length_simulation in length:
    
    savename = str(length_simulation)
    sim = SimM(length_simulation = length_simulation, advection_duration = advection_duration,
                              output_frequency = output_frequency, repeatdt = repeatdt, 
                              deltatime = deltatime, savename = savename,
                              domain = domain,
                              path_to_velocity = path_to_velocity,
                              path_to_grid = path_to_grid, data_in=data_in, bordercurrent = bordercurrent)
    sim.run_simulation()

#    beached_particles_list.append(sim.beached_particles)



beached_particles_list.append(sim.beached_particles)
length = [20,40,60,80,100,120,140,160,200,240,280,320]
plt.plot(length, beached_particles_list)
plt.ylabel('percentage particles on land')
plt.xlabel('simulation length (days)')

#%%


length_simulation = 50
bordercurrentlist = [0, 0.0000003, 0.0000006,0.0000009,0.0000012]
beached_particles_per_border_current = []
for bordercurrent in bordercurrentlist:
    
    savename = 'bc is,' + str(bordercurrent)
    sim = SimM(length_simulation = length_simulation, advection_duration = advection_duration,
                              output_frequency = output_frequency, repeatdt = repeatdt, 
                              deltatime = deltatime, savename = savename,
                              domain = domain,
                              path_to_velocity = path_to_velocity,
                              path_to_grid = path_to_grid, data_in=data_in, bordercurrent = bordercurrent)
    sim.run_simulation()
    Tm_obj_wob = TraM(sim, beaching_timescale)
    Tm_obj_wob.construct_matrix()
    beached_particles_per_border_current.append(sim.beached_particles)
    

plt.figure()
plt.plot(bordercurrentlist, beached_particles_per_border_current)
plt.axvline(x=0.0000012, color = 'red', label='line at x = {}'.format(0.0000012))
plt.ylabel('percentage particles on land')
plt.xlabel('value border current')
plt.legend()

#%%

Tm_obj_wob.construct_matrix()

#%%
#############plot trajectories####################

TP(sim_wbc).plot_trajectories()
TP(sim_wobc).plot_trajectories()





#%%
Tm_obj_wb.construct_matrix()
Tm_obj_wob.construct_matrix()

tm_wb = Tm_obj_wb.transition_matrix
tm_wob = Tm_obj_wob.transition_matrix

plt.figure()
plt.imshow(tm_wb)

plt.figure()
plt.imshow(tm_wob)

print('with bc:',sim_wbc.beached_particles)
print('without bc:', sim_wobc.beached_particles)
#%%

#Sim_data = SimulationObject.trajectory_data
#disttoshore = Sim_data['distancetoshore'].data


Sim_data = sim_wbc.trajectory_data
land = Sim_data['beached'].data
land_num = np.nan_to_num(land)
AllParticles = np.sum(land_num, axis = 1)


BeachedBarticles = np.count_nonzero(AllParticles)
#BeachedPercentage = BeachedBarticles/len(AllParticles)
print(BeachedPercentag






















#Sim_data = SimulationObject.trajectory_data
#disttoshore = Sim_data['distancetoshore'].data
land = Sim_data['beached'].data
land = np.nan_to_num(land)
AllParticles = np.sum(land, axis = 1)
BeachedBarticles = np.count_nonzero(AllParticles)
BeachedPercentage = BeachedBarticles/len(AllParticles)
print(BeachedPercentage)
