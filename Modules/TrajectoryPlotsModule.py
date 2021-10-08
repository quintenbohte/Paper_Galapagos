# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 15:07:51 2021

@author: quint
"""

import numpy as np
import xarray as xr
from SimulationModule import Simulation as SM
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

class TrajectoryPlots:
    
    
    def __init__(self, simulation_object = None, savename = None, domain = None):
        
        if simulation_object:
            self.simulation_object = simulation_object
            self.domain = self.simulation_object.domain
      
        else:
            self.domain = domain
        
        self.savename = savename
        self.lonCentreGrid = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\lonCentreGrid.npy')
        self.latCentreGrid = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\latCentreGrid.npy')
        self.lonCornerGrid = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\lonCornerGrid.npy')
        self.latCornerGrid = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\latCornerGrid.npy')
        self.bordercurrent = simulation_object.bordercurrent
        
        self.lonCentre = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\lonCentre.npy')
        self.latCentre = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\latCentre.npy')
        self.lonCorner = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\lonCorner.npy')
        self.latCorner = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\latCorner.npy')

        
        self.uvel = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\uvel.npy')
        self.vvel = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\vvel.npy')

        self.uvel_mask = np.where(np.isnan(self.uvel), 0, self.uvel)[0,:,:]
        self.vvel_mask = np.where(np.isnan(self.vvel), 0, self.vvel)[0,:,:]

        self.x_corners_grid = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\x_corners_grid.npy')
        self.y_corners_grid = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\y_corners_grid.npy')

        self.dx = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\dx.npy')
        self.dy = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\dy.npy')

        self.U_xloc = None
        self.V_xloc = None
        self.U_yloc = None
        self.V_yloc = None
    
        self.landmaskGalapagos = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\Landmask.npy')
        self.coastgrids = np.roll(np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\Coastgrids.npy'),1,0)
        self.x_plot = None
        self.y_plot = None
        
    def getclosest_ij(self,lats,lons,latpt,lonpt):    

        dist_lat = (lats-latpt)**2          # find squared distance of every point on grid
        dist_lon = (lons-lonpt)**2
        minindex_lat = dist_lat.argmin()    # 1D index of minimum dist_sq element
        minindex_lon = dist_lon.argmin()
        
        return minindex_lat, minindex_lon   # Get 2D index for latvals and lonvals arrays from 1D index

    
    def construct_v_and_u_locations(self):
        

        iy_min, ix_min = self.getclosest_ij(self.latCentre, self.lonCentre, self.domain[2], self.domain[0])
        iy_max, ix_max = self.getclosest_ij(self.latCentre, self.lonCentre, self.domain[3], self.domain[1])

        self.U_xloc = self.lonCornerGrid[iy_min:iy_max,ix_min:ix_max]
        self.U_yloc = self.latCentreGrid[iy_min:iy_max,ix_min:ix_max]
        
        self.V_xloc = self.lonCentreGrid[iy_min:iy_max,ix_min:ix_max]
        self.V_yloc = self.latCornerGrid[iy_min:iy_max,ix_min:ix_max]


    
    def load_trajectory_data(self):
        
        if self.savename:
            plot_data = xr.open_dataset(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Simualtion_' + self.savename + '.nc')
            self.x_plot = plot_data['lon'].T
            self.y_plot = plot_data['lat'].T
            
        else:
            plot_data = self.simulation_object.trajectory_data
            self.x_plot = plot_data['lon'].T
            self.y_plot = plot_data['lat'].T
            
                        
    def plot_trajectories(self):
        
        self.load_trajectory_data()
        
        cmap = plt.get_cmap('Blues')
        my_cmap = cmap(np.arange(cmap.N))
        my_cmap[:,-1] = 0
        my_cmap = ListedColormap(my_cmap)
        
        'Plotting the c-grid with the velocity vectors'
        plt.figure(figsize = (15,15))
        plt.axes()
        plt.title('', fontsize = 15)
        plt.xlim(-92,-88)
        plt.ylim(-2, 2)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('Longitude ($^\circ$)', fontsize = 25)
        plt.ylabel('Latitude ($^\circ$)', fontsize = 25)
        plt.title('bordercurrent: ' + str(self.bordercurrent), fontsize = 25)
        
        #landmask
        plt.pcolormesh(self.x_corners_grid, self.y_corners_grid, self.landmaskGalapagos,cmap='Blues',edgecolors='k',linewidth=1) 
        plt.pcolormesh(self.x_corners_grid, self.y_corners_grid, self.landmaskGalapagos,cmap='GnBu',edgecolors='k',linewidth=0.5) 
    
        #trajectroy
        plt.plot(self.x_plot, self.y_plot, linewidth = 1, marker = 'o', markersize = 3) #trajectories op particles
        
                
    def plot_velocity_vectors(self):
                
        self.construct_v_and_u_locations()
        
        cmap = plt.get_cmap('Blues')
        my_cmap = cmap(np.arange(cmap.N))
        my_cmap[:,-1] = 0
        my_cmap = ListedColormap(my_cmap)
        
        'Plotting the c-grid with the velocity vectors'
        plt.figure(figsize = (15,15))
        plt.axes()
        plt.title('', fontsize = 15)
        plt.xlim(-92,-88)
        plt.ylim(-2, 2)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('Longitude ($^\circ$)', fontsize = 25)
        plt.ylabel('Latitude ($^\circ$)', fontsize = 25)
        plt.title('Release Locations Particles', fontsize = 35)
        
        #landmask
        plt.pcolormesh(self.x_corners_grid, self.y_corners_grid, self.landmaskGalapagos,cmap='Blues',edgecolors='k',linewidth=1) 
        plt.pcolormesh(self.x_corners_grid, self.y_corners_grid, self.landmaskGalapagos,cmap='GnBu',edgecolors='k',linewidth=0.5) 
    
        #Velocity Field
        plt.scatter(self.U_xloc, self.U_yloc, c = self.uvel_mask,cmap='seismic',vmin=-0.1,vmax=0.1,edgecolor='k', label = 'U') #U vectors
        plt.scatter(self.V_xloc, self.V_yloc, c = self.vvel_mask, cmap='PRGn',vmin=-0.1,vmax=0.1,edgecolor='k', label = 'V') #V vectors
        plt.quiver(self.U_xloc,self.U_yloc,self.uvel[0,:,:],np.zeros(self.U_xloc.shape),angles='xy', scale_units='xy', scale=6, width=0.005) #vector arrows
        plt.quiver(self.V_xloc,self.V_yloc, np.zeros(self.V_xloc.shape), self.vvel[0,:,:], angles='xy', scale_units='xy', scale=6, width=0.005) #vector arrows
    
                  
        
#%%
         
#length_simulation = 10 #unit: days (for how long do we deploy particles)
#advection_duration = 60 #unit: days (how long does one particle advect in the fields)
#output_frequency = 1     #unit: hours
#repeatdt = 48         #unit: hours
#deltatime = 1           #dt in hours
#bordercurrent = 0.00000000
#domain = [-92, -88, -2, 2]
#
#data_in = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input'
#savename = 'withoutBordercurrent'
#path_to_velocity = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input\RGEMS_2010_Surf.nc'
#path_to_grid = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input\RGEMS3_Surf_grid.nc'
#
#
#sim_wobc = SM(length_simulation = length_simulation, advection_duration = advection_duration,
#                              output_frequency = output_frequency, repeatdt = repeatdt, 
#                              deltatime = deltatime, savename = savename,
#                              domain = domain,
#                              path_to_velocity = path_to_velocity,
#                              path_to_grid = path_to_grid, data_in=data_in, bordercurrent = bordercurrent)
#
#
#sim_wobc.run_simulation()  
##%%
#      
#
#savename = 'withoutBordercurrent'
#            
#TrajectoryPlots(sim_wobc).plot_trajectories()
#TrajectoryPlots(sim_wobc).plot_velocity_vectors()
#
#






    
    