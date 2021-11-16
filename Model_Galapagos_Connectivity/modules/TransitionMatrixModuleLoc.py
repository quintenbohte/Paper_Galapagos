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
        
        self.parent_dir = str(Path(os.path.abspath(__file__)).parents[1])
        self.gridsdataframe = pd.read_csv(self.parent_dir + '\data\input\galapagos_field_data\GridsDataFrame.csv')
        self.beaching_timescale = beaching_timescale  
        self.coastcells = traj_data_dict.get('coastgrids')
        self.gridnumbers = traj_data_dict.get('gridnumber')
        self.dt = 60*delta_time
        self.savename = savename
        
        self.trajectory_lons = None
        self.trajectory_lats = None
        self.transition_matrix =  None
    
#    
#    def TrajectroyOutputLocations(self):
#    
#        trajectory_lons = self.trajectory_data.variables['lon'][:]
#        trajectory_lats = self.trajectory_data.variables['lat'][:]
#        
#        self.trajectory_lons = np.nan_to_num(np.asarray(trajectory_lons), nan = 1000)
#        self.trajectory_lats = np.nan_to_num(np.asarray(trajectory_lats), nan = 1000)

    def BeachProb(self):
  
        Prob = math.exp((-self.dt*60)/(self.beaching_timescale*24*60*60))
        
        return Prob
    
    
    def CreateEmptyTransMatrix(self):
    
        self.transition_matrix = np.zeros((len(self.gridsdataframe), len(self.gridsdataframe)))
#        self.transition_matrix = np.fill_diagonal( self.transition_matrix, np.nan)
        
    def save_transition_matrix(self):
        
#        np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\TransitionMatrices\TM' + self.savename ,self.transition_matrix)
#        np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Model_local_pc\data\output\transition_matrices\Tm' + self.savename , self.transition_matrix)
        np.save(self.parent_dir + '\data\output\Transition_matrices\Tm' + self.savename , self.transition_matrix)
        
    def construct_matrix(self):
        
        count = 0
#        self.TrajectroyOutputLocations()
        self.CreateEmptyTransMatrix()
        
        Timesteparray = np.zeros(len(self.coastcells[0,:]))
        
        for particle in range(self.coastcells.shape[0]):
#            print(particle)
            breaker = False
            print(self.coastcells.shape[0]-particle, 'iterations to go',)
            
            for Timestep in range(len(self.coastcells[particle,:])):
#                print(Timestep)
                CoastcellValue = self.coastcells[particle,Timestep]
                
                if CoastcellValue == 1:
                    Prob = self.BeachProb()
                    
                    if random.uniform(0,1) > Prob:
#        
#                        print('particle:', particle)
#                        print('timestep:', Timestep)
#                        print('gridnumber:', self.gridnumbers[particle,Timestep])
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
                    
        self.save_transition_matrix()
        
        print(count, 'particles have beached, which is', np.round(count/self.coastcells.shape[0], 2)*100, 'percent')    

        


















