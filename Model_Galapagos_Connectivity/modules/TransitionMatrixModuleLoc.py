# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 21:06:24 2021

@author: quint
"""


import numpy as np
import pandas as pd
import math 
from parcels import rng as random
from SimulationModule import Simulation as SM
import os
from pathlib import Path

class ConstructTransitionMatrix:
    
    def __init__(self, simulation_dict, beaching_timescale):
        
        self.dict_simulation = simulation_dict
        self.coastcells = simulation_dict['trajectory_data']['coastcell'].data
        self.gridnumbers = simulation_dict['trajectory_data']['gridnumber'].data
        self.gridsdataframe = pd.read_csv(r'C:\Users\quint\Documents\Quinten_studie\Master_Thesis\Data\Output_Data\ReleaseLocations\GridsDataFrame.csv')
        self.trajectory_data = simulation_dict['trajectory_data']
        self.dt = 60*simulation_dict['deltatime']
        self.beaching_timescale = beaching_timescale  
        self.savename = simulation_dict['savename']
        
        self.trajectory_lons = None
        self.trajectory_lats = None
        self.transition_matrix =  None
    
    
    def TrajectroyOutputLocations(self):
    
        trajectory_lons = self.trajectory_data.variables['lon'][:]
        trajectory_lats = self.trajectory_data.variables['lat'][:]
        
        self.trajectory_lons = np.nan_to_num(np.asarray(trajectory_lons), nan = 1000)
        self.trajectory_lats = np.nan_to_num(np.asarray(trajectory_lats), nan = 1000)

    def BeachProb(self):
  
        Prob = math.exp((-self.dt*60)/(self.beaching_timescale*24*60*60))
        
        return Prob
    
    
    def CreateEmptyTransMatrix(self):
    
        self.transition_matrix = np.zeros((len(self.gridsdataframe), len(self.gridsdataframe)))
#        self.transition_matrix = np.fill_diagonal( self.transition_matrix, np.nan)
        
    def save_transition_matrix(self):
        
        np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\TransitionMatrices\TM' + self.savename ,self.transition_matrix)

    def construct_matrix(self):
        
        count = 0
        self.TrajectroyOutputLocations()
        self.CreateEmptyTransMatrix()
        
        Timesteparray = np.zeros(len(self.coastcells[0,:]))
        
        for particle in range(self.trajectory_lons.shape[0]):
        
            breaker = False
            print(self.trajectory_lons.shape[0]-particle, 'iterations to go',)
            
            for Timestep in range(len(self.coastcells[particle,:])):
        #        print(Timestep)
                CoastcellValue = self.coastcells[particle,Timestep]
        
                if CoastcellValue == 1:
                    Prob = self.BeachProb()
        
                    if random.uniform(0,1) > Prob:
        
                        EndGrid = int(self.gridnumbers[particle, Timestep])
        
                        if particle < 393:
                            BeginGrid = int(particle)
                            
                        else:
                            BeginGrid = int(particle%393)
                        
                        print('begin grid is', BeginGrid)
                        print('end grid is', EndGrid)
                        print('   ')
        
                        count += 1
                        
        
                        self.transition_matrix[BeginGrid,EndGrid-1] += 1
                        
                        Timesteparray[Timestep] += 1
                        breaker = True
                        
                        break
        print(count, 'particles have beached, which is', np.round(count/self.trajectory_lons.shape[0], 2)*100, 'percent')    
        self.save_transition_matrix()
        
class ConstructTransitionMatrixOnlyTrajData:
    
    
    def __init__(self, traj_data_dict, delta_time, savename, beaching_timescale):
        
        #parent directory
        self.parent_dir = str(Path(os.path.abspath(__file__)).parents[1])
        
        #gridsdataframe, this dataframe includes all information about the coastgrids
        self.gridsdataframe = pd.read_csv(self.parent_dir + '/data/input/galapagos_field_data/GridsDataFrame.csv')
        self.beaching_timescale = beaching_timescale  
        
        #the coastcells dataframe, is a matrix with the same shape as the particle trajectories. It contains all timesteps for each particle and for each timestep
        #it gives a 1 when the particle is in a coastgrid (potential beaching location) and zero otherwise
        self.coastcells = traj_data_dict.get('coastgrids')
        
        #the gridnumber matrix, is a matrix with each row representing a particle and each columns representing a timestep. The columns values are the gridnumbers of the coastgrids.
        #So this matrix tells us in what coastgrid the particle is a what timestep
        self.gridnumbers = traj_data_dict.get('gridnumber')
        self.dt = 60*delta_time
        self.savename = savename
        
        self.trajectory_lons = None
        self.trajectory_lats = None
        self.transition_matrix =  None
    

    #function that computes a beaching probability based on the beaching timescale.
    def BeachProb(self):
  
        Prob = math.exp((-self.dt*60)/(self.beaching_timescale*24*60*60))
        
        return Prob
    
    #function that makes the empty transition matrix, which will be filled with values. The row numbers represnt the release grid and the columns numbers
    #represent the beach grids
    def CreateEmptyTransMatrix(self):
    
        self.transition_matrix = np.zeros((len(self.gridsdataframe), len(self.gridsdataframe)))

    #function that saves the transition matrix
    def save_transition_matrix(self):
        
        np.save(self.parent_dir + '/data/output/Transition_matrices/Tm' + self.savename , self.transition_matrix)
        
    #function that construct the transition matrux
    def construct_matrix(self):
        
        count = 0
        
        #create empty transition matrix
        self.CreateEmptyTransMatrix()
        
        #construct the timestep array. Each entry of this array represents a timestep of the simulation. Everytime a particle beaches, we add 1 to the entry of the timestep that
        #corresponds to the timestep when the particle. In this manner we know how many particles have beached on what timestep, this may be usefull to information when we want to
        #investigate what influence the beaching timescale has on teh duration a particle beaches. It may also give us an idea in what time range most particles beaches
        Timesteparray = np.zeros(len(self.coastcells[0,:]))
        
        #loop over all different particles
        for particle in range(self.coastcells.shape[0]):

            breaker = False
            print(self.coastcells.shape[0]-particle, 'iterations to go',)
            
            #loop over all timesteps of a particle.
            for Timestep in range(len(self.coastcells[particle,:])):

                #check the coastcell value. 
                CoastcellValue = self.coastcells[particle,Timestep]
                
                #if the coastcell value of the particle is 1, we know that it is in a coastgrid and we should compute the beaching probability.
                if CoastcellValue == 1:
                    #compute the beaching probability
                    Prob = self.BeachProb()
                    
                    #compute a random number between 0-1, when this random number is larger than the beaching probability, we mark the particle as beached.
                    if random.uniform(0,1) > Prob:
#        
                        #when the particle is marked as beached, define the end grid. This end grid is the gridnumber of the coastgrid where it beached. 
                        EndGrid = int(self.gridnumbers[particle, Timestep])
                        
                        
                        #since we have 393 release locations, the gridnumber of the begingrid is equal to the particle number when particle number is lower than 393
                        if particle < 393:
                            BeginGrid = int(particle)
                            
                        #when the particle number is larger than 393, we use the % operator and divide it by 393, the result is the begin grid.
                        else:
                            BeginGrid = int(particle%393)
                        
                        print('begin grid is', BeginGrid)
                        print('end grid is', EndGrid)
                        print('   ')
        
                        count += 1
                        
                        #we now determined the release and beach grid. We add 1 to the entry corresponding to the begin and end grid. 
                        self.transition_matrix[BeginGrid,EndGrid-1] += 1
                        
                        #we add 1
                        Timesteparray[Timestep] += 1
                        breaker = True
                        
                        break
                    
#        self.save_transition_matrix()
        
        print(count, 'particles have beached, which is', np.round(count/self.coastcells.shape[0], 2)*100, 'percent')    

        
#
#class ConstructTransitionMatrixOnlyTrajData:
#    
#    
#    
#    def __init__(self, traj_data_dict, delta_time, savename, beaching_timescale):
#        
#        self.parent_dir = str(Path(os.path.abspath(__file__)).parents[1])
#        self.gridsdataframe = pd.read_csv(self.parent_dir + '\data\input\galapagos_field_data\GridsDataFrame.csv')
#        self.beaching_timescale = beaching_timescale  
#        self.coastcells = traj_data_dict.get('coastgrids')
#        self.gridnumbers = traj_data_dict.get('gridnumber')
#        self.dt = 60*delta_time
#        self.savename = savename
#        
#        self.trajectory_lons = None
#        self.trajectory_lats = None
#        self.transition_matrix =  None
#    
##    
##    def TrajectroyOutputLocations(self):
##    
##        trajectory_lons = self.trajectory_data.variables['lon'][:]
##        trajectory_lats = self.trajectory_data.variables['lat'][:]
##        
##        self.trajectory_lons = np.nan_to_num(np.asarray(trajectory_lons), nan = 1000)
##        self.trajectory_lats = np.nan_to_num(np.asarray(trajectory_lats), nan = 1000)
#
#    def BeachProb(self):
#  
#        Prob = math.exp((-self.dt*60)/(self.beaching_timescale*24*60*60))
#        
#        return Prob
#    
#    
#    def CreateEmptyTransMatrix(self):
#    
#        self.transition_matrix = np.zeros((len(self.gridsdataframe), len(self.gridsdataframe)))
##        self.transition_matrix = np.fill_diagonal( self.transition_matrix, np.nan)
#        
#    def save_transition_matrix(self):
#        
##        np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\TransitionMatrices\TM' + self.savename ,self.transition_matrix)
##        np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Model_local_pc\data\output\transition_matrices\Tm' + self.savename , self.transition_matrix)
#        np.save(self.parent_dir + '\data\output\Transition_matrices\Tm' + self.savename , self.transition_matrix)
#        
#    def construct_matrix(self):
#        
#        count = 0
##        self.TrajectroyOutputLocations()
#        self.CreateEmptyTransMatrix()
#        
#        Timesteparray = np.zeros(len(self.coastcells[0,:]))
#        
#        for particle in range(self.coastcells.shape[0]):
##            print(particle)
#            breaker = False
#            print(self.coastcells.shape[0]-particle, 'iterations to go',)
#            
#            for Timestep in range(len(self.coastcells[particle,:])):
##                print(Timestep)
#                CoastcellValue = self.coastcells[particle,Timestep]
#                
#                if CoastcellValue == 1:
#                    Prob = self.BeachProb()
#                    
#                    if random.uniform(0,1) > Prob:
##        
##                        print('particle:', particle)
##                        print('timestep:', Timestep)
##                        print('gridnumber:', self.gridnumbers[particle,Timestep])
#                        EndGrid = int(self.gridnumbers[particle, Timestep])
#                        
#                        
#                        
#                        if particle < 393:
#                            BeginGrid = int(particle)
#                            
#                        else:
#                            BeginGrid = int(particle%393)
#                        
#                        print('begin grid is', BeginGrid)
#                        print('end grid is', EndGrid)
#                        print('   ')
#        
#                        count += 1
#                        
#        
#                        self.transition_matrix[BeginGrid,EndGrid-1] += 1
#                        
#                        Timesteparray[Timestep] += 1
#                        breaker = True
#                        
#                        break
#                    
##        self.save_transition_matrix()
#        
#        print(count, 'particles have beached, which is', np.round(count/self.coastcells.shape[0], 2)*100, 'percent')    
#
#        


















