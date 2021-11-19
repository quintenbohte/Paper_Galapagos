# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 10:38:06 2021

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



#from SimpleFunctions import Functions as Func

############ FUNCTIONS ####################

def Convert_to_single_particles(data, length_simulation, advection_duration, output_frequency, repeatdt, deltatime):

    advection_duration_hours = int(advection_duration*(24/output_frequency))+1 
    advected_timesteps = len(data[1,:])
    number_particles = math.ceil(advected_timesteps/advection_duration_hours)
    ReleaseLon = list(np.load(parent_dir + '\data\input\galapagos_field_data\ReleaseLon.npy'))

    release_locations = len(ReleaseLon)

    output_trajectory = np.zeros(((len(data[:,1])*number_particles) ,advection_duration_hours))

    for trajectory_number in range(len(data)):
        
        trajectory = data[trajectory_number,:]
    
        remaining_time_steps = len(trajectory)%advection_duration_hours
    
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
                
                output_trajectory[index,0:len(last_particle)] = last_particle
                
    return output_trajectory

############ RELEASE LOCATIONS ####################


ReleaseLat = list(np.load(parent_dir + '\data\input\galapagos_field_data\ReleaseLat.npy'))
ReleaseLon = list(np.load(parent_dir + '\data\input\galapagos_field_data\ReleaseLon.npy'))


############ SIMULATIONS SPECIFICATIONS ####################

length_simulation = 100 #unit: days (for how long do we deploy particles)
advection_duration = 60 #unit: days (how long does one particle advect in the fields)
output_frequency = 1     #unit: hours
repeatdt = 24        #unit: hours
deltatime = 1           #dt in hours

#data_in = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input'
#data_in = '\model_local_pc\da'
savename = 'test'

domain = [-92, -88, -2, 2]

############ GET INDICES ####################
def getclosest_ij(lats,lons,latpt,lonpt):    
    """Function to find the index of the closest point to a certain lon/lat value."""
    dist_lat = (lats-latpt)**2          # find squared distance of every point on grid
    dist_lon = (lons-lonpt)**2
    minindex_lat = dist_lat.argmin()    # 1D index of minimum dist_sq element
    minindex_lon = dist_lon.argmin()
    return minindex_lat, minindex_lon   # Get 2D index for latvals and lonvals arrays from 1D index


dfile = Dataset(parent_dir + '\data\input\grid_data\RGEMS3_Surf_grid.nc')
lon = dfile.variables['XC'][:]
lat = dfile.variables['YC'][:]
loncor = dfile.variables['XG'][:]
latcor = dfile.variables['YG'][:]
iy_min, ix_min = getclosest_ij(lat, lon, domain[2], domain[0])
iy_max, ix_max = getclosest_ij(lat, lon, domain[3], domain[1])

londomain = lon[ix_min:ix_max]
latdomain = lat[iy_min:iy_max]

loncordomain = lon[ix_min:ix_max]
latcordomain = lat[iy_min:iy_max]
############ CREATE FIELDSET ####################

#PathToVelocity = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input\velocity_data_2008_2009\*'
PathToVelocity = parent_dir + '\data\input\RGEMS_surf\RGEMS_surf_all_years\*'
#PathToVelocityrel = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Model_local_pc\data\input\velocity_data\velocity_data_all_years\*'
#PathToGrid = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input\RGEMS3_Surf_grid.nc'
PathToGrid = parent_dir + '\data\input\grid_data\RGEMS3_Surf_grid.nc'
#PathToGrid = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Model_local_pc\data\input\grid_data\RGEMS3_Surf_grid.nc'


VelocityField = sorted(glob(PathToVelocity))
GridField = glob(PathToGrid)


files_MITgcm = {'U': {'lon': GridField, 'lat': GridField, 'data': VelocityField},                                             
            'V': {'lon': GridField, 'lat': GridField, 'data': VelocityField}}     

variables = {'U': 'UVEL', 'V': 'VVEL'}
dimensions = {'lon': 'XG', 'lat': 'YG', 'time': 'time'}

indices = {'lon': range(ix_min,ix_max), 
           'lat': range(iy_min,iy_max)}


fieldset = FieldSet.from_mitgcm(files_MITgcm,
                                     variables,
                                     dimensions,
                                     indices = indices)


############ ADD EXTRA FIELDS/CONSTANTS TO FIELDSET ####################

#load fields

landmask = np.load(parent_dir + '\data\input\galapagos_field_data\Landmask.npy')
coastgrids = np.load(parent_dir + '\data\input\galapagos_field_data\Coastgrids.npy')
coastgrids = np.roll(coastgrids,1,0)
gridnumbermask = np.load(parent_dir + '\data\input\galapagos_field_data\GridNumberMask.npy')
gridnumbermask = np.roll(gridnumbermask,1,0)

fieldset.add_constant('advection_duration',advection_duration)    
fieldset.add_constant('lon_max',lon[ix_max] - 0.2)
fieldset.add_constant('lon_min',lon[ix_min] + 0.2)
fieldset.add_constant('lat_max',lat[iy_max] - 0.2)
fieldset.add_constant('lat_min',lat[iy_min] + 0.2)


fieldset.add_field(Field('landmask',
                         data = landmask,
                         lon = londomain,
                         lat = latdomain,
                         mesh='spherical',
                         interp_method = 'nearest'))

fieldset.add_field(Field('coastgrids',
                         data = coastgrids,
                         lon = londomain,
                         lat = latdomain,
                         mesh='spherical',
                         interp_method = 'nearest'))

fieldset.add_field(Field('gridnumbermask',
                         data = gridnumbermask,
                         lon = londomain,
                         lat = latdomain,
                         mesh='spherical',
                         interp_method = 'nearest'))



#BorU = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\borderCurrentU.npy')
#BorUlon = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\BorderUlon.npy')
#BorUlat = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\BorderUlat.npy')
#
#BorV = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\borderCurrentV.npy')
#BorVlon = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\BorderVlon.npy')
#BorVlat = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\BorderVlat.npy')

#fieldset.add_field(Field('BorderCurrentU', BorU, BorUlon, BorUlat))
#fieldset.add_field(Field('BorderCurrentV', BorV, BorVlon, BorVlat))


DistanceFromShore = np.load(parent_dir + '\data\input\galapagos_field_data\distance2shore.npy')
x = np.linspace(domain[0],domain[1],len(londomain))
y = np.linspace(domain[2],domain[3],len(londomain))
lon, lat = np.meshgrid(x,y)
fieldset.add_field(Field('distance2shore', DistanceFromShore, lon, lat))


############ ADD KERNELS ####################


def AdvectionRK4(particle, fieldset, time):
    """ Only advect particles that are not out of bounds"""
    
    if (particle.lon < fieldset.lon_max and
        particle.lon > fieldset.lon_min and
        particle.lat < fieldset.lat_max and
        particle.lat > fieldset.lat_min):
    
        (u1, v1) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
        lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
#        particle.u1p = u1
#        particle.v1p = v1
        if (lon1 < fieldset.lon_max and
            lon1 > fieldset.lon_min and
            lat1 < fieldset.lat_max and
            lat1 > fieldset.lat_min):

            (u2, v2) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat1, lon1]
            lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)

#            particle.u2p = u2
#            particle.v2p = v2
            if (lon2 < fieldset.lon_max and
                lon2 > fieldset.lon_min and
                lat2 < fieldset.lat_max and
                lat2 > fieldset.lat_min):

                (u3, v3) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat2, lon2]
                lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)

#                particle.u3p = u3
#                particle.v3p = v3

                if (lon3 < fieldset.lon_max and
                    lon3 > fieldset.lon_min and
                    lat3 < fieldset.lat_max and
                    lat3 > fieldset.lat_min):

                    (u4, v4) = fieldset.UV[time + particle.dt, particle.depth, lat3, lon3]
                    
#                    particle.u4p = u4
#                    particle.v4p = v4
#                    
#                    particle.displ_x = (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
#                    particle.displ_y = (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
                    
                
##                    if particle.displ_x == math.nan:
##                        
##                        particle.lon = particle.lon
##                        particle.lat = particle.lat
##                        
##                    else:
#                    if (particle.u1p == math.nan or 
#                       particle.u2p == math.nan or
#                       particle.u3p == math.nan or
#                       particle.u4p == math.nan or
#                       particle.v1p == math.nan or
#                       particle.v2p == math.nan or
#                       particle.v3p == math.nan or
#                       particle.v4p == math.nan):
#                        
#                        particle.lon = particle.lon
#                        particle.lat = particle.lat
#                        
#                        
#                        
#                        
#                    else:
#                    
                    particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
                    particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
                    
                    
#    else:
#        particle.lon = particle.lon
#        particle.lat = particle.lat


def Age(fieldset, particle, time):
    """ Delete particles when reaching age specified by advection_duration """
    
    particle.age = particle.age + particle.delta_time*60*60
    
  
#    if particle.age > fieldset.advection_duration*86400:
#        particle.delete()

def SampleInfo(fieldset, particle, time):
    if (particle.lon < fieldset.lon_max and
        particle.lon > fieldset.lon_min and
        particle.lat < fieldset.lat_max and
        particle.lat > fieldset.lat_min):
        
        particle.distance = fieldset.distance[time, particle.depth, particle.lat, particle.lon]
        particle.island = fieldset.island[time, particle.depth, particle.lat, particle.lon]    

#def DeleteParticle(particle, fieldset, time):
#    if particle.lon < fieldset.lon_min+0.1 or particle.lon > fieldset.lon_max-0.1 or particle.lat < fieldset.lat_min+0.1 or particle.lat > fieldset.lat_max-0.1: 
#        particle.delete()
    

def beachtesting(particle, fieldset, time):
    
    landcheck = fieldset.landmask[time, particle.depth, particle.lat, particle.lon]
    
    if landcheck == 1: 

        particle.beached = 1
        
def coasttesting(particle, fieldset, time):
    
    if (particle.lon < fieldset.lon_max and
        particle.lon > fieldset.lon_min and
        particle.lat < fieldset.lat_max and
        particle.lat > fieldset.lat_min):
    
        particle.coastcell = fieldset.coastgrids[time, particle.depth, particle.lat, particle.lon]
        
def gridnumbertesting(particle, fieldset, time):
    
    if (particle.lon < fieldset.lon_max and
        particle.lon > fieldset.lon_min and
        particle.lat < fieldset.lat_max and
        particle.lat > fieldset.lat_min):
    
        particle.gridnumber = fieldset.gridnumbermask[time, particle.depth, particle.lat, particle.lon]
        
        
def delete_particle(particle, fieldset, time):
   # This delete particle format from Philippe Delandmeter
   # https://github.com/OceanParcels/Parcelsv2.0PaperNorthSeaScripts/blob/master/northsea_mp_kernels.py
   print("Particle [%d] lost !! (%g %g %g %g)" % (
   particle.id, particle.lon, particle.lat, particle.depth, particle.time))
   particle.delete()
   
#def AntiBeachNudging(particle,fieldset,time):
#    """    
#    The nudging current is 1 m s^-1, which ought to be sufficient to overpower
#    any coastal current (I hope) and push our particle back out to sea so as to
#    not get stuck
#    
#    update 11/03/2020: Following tests and discussions with Cleo, the nudging 
#    current will now kick in starting at 500m from the coast, since otherwise 
#    the particles tended to get stuck if we used the velocity treshhold. 
#    """
#    particle.distancetoshore = fieldset.distance2shore[time,particle.depth,particle.lat,particle.lon]
#    if particle.distancetoshore < 2:
#        print('test')
#        borUab,borVab=fieldset.BorderCurrentU[time, particle.depth, particle.lat, particle.lon],fieldset.BorderCurrentV[time, particle.depth, particle.lat, particle.lon]
#        particle.lon += 0.0000000*borUab*particle.dt
#        particle.lat += 0.0000000*borVab*particle.dt

def SetParticleBack(particle,fieldset,time):
    
    if particle.age > fieldset.advection_duration*86400:
        
        print('test')
        particle.lon = particle.startlon
        particle.lat = particle.startlat
        particle.age = 0
        
        
        
        
############ CREATE PARTICLE ####################

class GalapagosParticle(JITParticle):
    distance = Variable('distance', initial = 0) 
    age = Variable('age', dtype=np.float32, initial = 0.)
    beached = Variable('beached', dtype=np.int32, initial = 0.) 
    coastcell = Variable('coastcell', dtype=np.int32, initial = 0.)
    gridnumber = Variable('gridnumber', dtype=np.int32, initial = 0.)
    distancetoshore = Variable('distancetoshore', dtype=np.int32, initial = 0.)
    startlon = Variable('startlon', dtype=np.float32, initial = ReleaseLon)
    startlat = Variable('startlat', dtype=np.float32, initial = ReleaseLat)
    delta_time = Variable('delta_time', dtype=np.float32, initial = deltatime)
    
#    u4_test =  Variable('u4_test', dtype=np.float32, initial = 0.)
#    displ_x = Variable('displ_x', dtype=np.float32, initial = 0.)
#    displ_y = Variable('displ_y', dtype=np.float32, initial = 0.)
#    u1p =  Variable('u1p', dtype=np.float32, initial = 0.)
#    u2p =  Variable('u2p', dtype=np.float32, initial = 0.)
#    u3p =  Variable('u3p', dtype=np.float32, initial = 0.)
#    u4p =  Variable('u4p', dtype=np.float32, initial = 0.)
#    v1p =  Variable('v1p', dtype=np.float32, initial = 0.)
#    v2p =  Variable('v2p', dtype=np.float32, initial = 0.)
#    v3p =  Variable('v3p', dtype=np.float32, initial = 0.)
#    v4p =  Variable('v4p', dtype=np.float32, initial = 0.)

if repeatdt == None:
    print('NoRepeat')
    pset = ParticleSet(fieldset = fieldset,
                           pclass = GalapagosParticle,
                           lon = ReleaseLon,
                           lat = ReleaseLat,
                           time = 0)
else:
    print('repeat')
    pset = ParticleSet(fieldset = fieldset,
                           pclass = GalapagosParticle,
                           lon = ReleaseLon,
                           lat = ReleaseLat,
                           repeatdt = delta(hours = repeatdt),
                           time = 0)  

############ EXECUTION PARTICLE SET ####################


outfile = pset.ParticleFile(parent_dir + '\data\output\simulations\Simulation_' + savename + '.nc', outputdt = delta(hours = output_frequency))

kernels = (pset.Kernel(AdvectionRK4) + pset.Kernel(Age) + 
           pset.Kernel(beachtesting) + pset.Kernel(coasttesting) +
           pset.Kernel(gridnumbertesting)+ pset.Kernel(SetParticleBack))


pset.execute(kernels,
             runtime=delta(days=advection_duration),
             dt=delta(hours=deltatime),
             recovery = {ErrorCode.ErrorInterpolation : delete_particle},
             output_file=outfile)

pset.repeatdt = None

pset.execute(kernels,
             runtime=delta(days=length_simulation - advection_duration),
             dt=delta(hours=deltatime),
             recovery = {ErrorCode.ErrorInterpolation : delete_particle},
             output_file=outfile)


outfile.export()
outfile.close()

#Sim_data = xr.open_dataset(parent_dir + '\data\output\simulations\Simulation_' + savename + '.nc')
#lon_traj = Sim_data['lon'].data
#lat_traj = Sim_data['lat'].data 
#age_traj = Sim_data['age'].data 
#gridnumber_traj = Sim_data['gridnumber'].data
#coastgrid_traj = Sim_data['coastcell'].data



############################################## CONVERT DATA TO SINGLE PARTICLES ##############################################################

#We set the particles back after the advection time. Due to this, one row contains all seperate trajectrories of a particle. 
#Therfore we have to chop this rows into the single trajectories and convert this to a matrix with all rows representing a single particle
#trajectrory. This is done with the convert_to_single_particle function.


Sim_data = xr.open_dataset(parent_dir + '\data\output\simulations\Simulation_' + savename + '.nc')
lon_traj = Sim_data['lon'].data
lat_traj = Sim_data['lat'].data 
age_traj = Sim_data['age'].data 
gridnumber_traj = Sim_data['gridnumber'].data 
coastgrid_traj = Sim_data['coastcell'].data 

#convert all simulation data to single particles.
age_per_particle = Convert_to_single_particles(age_traj, length_simulation, advection_duration, output_frequency, repeatdt, deltatime)
lon_per_particle = Convert_to_single_particles(lon_traj, length_simulation, advection_duration, output_frequency, repeatdt, deltatime)
lat_per_particle = Convert_to_single_particles(lat_traj, length_simulation, advection_duration, output_frequency, repeatdt, deltatime)
grid_numb_per_particle = Convert_to_single_particles(gridnumber_traj, length_simulation, advection_duration, output_frequency, repeatdt, deltatime)
coastgrid_per_particle = Convert_to_single_particles(coastgrid_traj, length_simulation, advection_duration, output_frequency, repeatdt, deltatime)

#create the dictionary with all the simulation data in it. This dictionary is used to construct the transition matrices
trajectory_data = {
                    'age':age_per_particle,
                    'lon':lon_per_particle,
                   'lat':lat_per_particle,
                   'gridnumber':grid_numb_per_particle,
                   'coastgrids':coastgrid_per_particle}

#save the dictionary. 
with open(parent_dir + '\data\output\simulations\Simulation_' + savename +  '.dictionary', 'wb') as config_dictionary_file:
        
         pickle.dump(trajectory_data, config_dictionary_file)










#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#data = age_traj
#
#advection_duration_hours = int(advection_duration*(24/output_frequency))+1 
#advected_timesteps = len(data[1,:])
#number_particles = math.ceil(advected_timesteps/advection_duration_hours)
#ReleaseLon = list(np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\ReleaseLocations\ReleaseLon.npy'))
#
#release_locations = len(ReleaseLon)
#
#output_trajectory = np.zeros(((len(data[:,1])*(number_particles+100)), advection_duration_hours))
#
#for trajectory_number in range(len(data)):
#    
#    trajectory = data[trajectory_number,:]
#
#    remaining_time_steps = len(trajectory)%advection_duration_hours
#
#    full_particle_trajectories = trajectory[0:len(trajectory)-remaining_time_steps]
#    
#    seperate_particles = np.split(full_particle_trajectories,number_particles-1)
#    last_particle = trajectory[len(full_particle_trajectories):(len(full_particle_trajectories)+remaining_time_steps)]
#
#    particle_set = int(math.floor(trajectory_number/release_locations))
#
#    for particle_number in range((number_particles)): 
#        
#        if particle_number < number_particles-1:
#        
#            print('part_set:',particle_set)
#            indices_current_particle_set = particle_set * (release_locations*number_particles)
#            print('current_ind:',indices_current_particle_set)
#            index = indices_current_particle_set + ((trajectory_number%release_locations) + (release_locations*particle_number))
#            print('index:',index)
#            
#            output_trajectory[index,:] = seperate_particles[particle_number]
#            
#        else:
#            print('part_set:',particle_set)
#            indices_current_particle_set = particle_set * (release_locations*number_particles)
#            print('current_ind:',indices_current_particle_set)
#            index = indices_current_particle_set + ((trajectory_number%release_locations) + (release_locations*particle_number))
#            print('index:',index)
#            
#            output_trajectory[index,0:len(last_particle)] = last_particle
#
#
#age_per_particle = Convert_to_single_particles(age_traj, length_simulation, advection_duration, output_frequency, repeatdt, deltatime)
#lon_per_particle = Convert_to_single_particles(lon_traj, length_simulation, advection_duration, output_frequency, repeatdt, deltatime)
#lat_per_particle = Convert_to_single_particles(lat_traj, length_simulation, advection_duration, output_frequency, repeatdt, deltatime)
#grid_numb_per_particle = Convert_to_single_particles(gridnumber_traj, length_simulation, advection_duration, output_frequency, repeatdt, deltatime)
#coastgrid_per_particle = Convert_to_single_particles(coastgrid_traj, length_simulation, advection_duration, output_frequency, repeatdt, deltatime)
#
#
##%%
#trajectory_data = {'age':age_per_particle,
#                   'lon':lon_per_particle,
#                   'lat':lat_per_particle,
#                   'gridnumber':grid_numb_per_particle,
#                   'coastgrids':coastgrid_per_particle}
#
#with open(parent_dir + '\data\output\simulations\Simulation_' + savename +  '.dictionary', 'wb') as config_dictionary_file:
#        
#         pickle.dump(trajectory_data, config_dictionary_file)
#
#
#
#
#
#
#
#
#
##age_subset = age_traj[0:10000]
##
###age_per_particle = Convert_to_single_particles(age_traj, length_simulation, advection_duration, output_frequency, repeatdt, deltatime)
###lon_per_particle = Convert_to_single_particles(lon_traj, length_simulation, advection_duration, output_frequency, repeatdt, deltatime)
###lat_per_particle = Convert_to_single_particles(lat_traj, length_simulation, advection_duration, output_frequency, repeatdt, deltatime)
###grid_numb_per_particle = Convert_to_single_particles(gridnumber_traj, length_simulation, advection_duration, output_frequency, repeatdt, deltatime)
###coastgrid_per_particle = Convert_to_single_particles(coastgrid_traj, length_simulation, advection_duration, output_frequency, repeatdt, deltatime)
##
###age_subset = age_per_particle[0:10000]
##
###trajectory_data = {'age':age_per_particle,
###                   'lon':lon_per_particle,
###                   'lat':lat_per_particle,
###                   'gridnumber':grid_numb_per_particle,
###                   'coastgrids':coastgrid_per_particle}
###
###
###
###def test_convert_to_single_particle(matrix):
###    
###    for i in range(len(matrix[:,1])): #row
###        startage = matrix[i,0]
###        if startage != 0:
###            print('start age')
###            print(i)
###            
###        for j in range(len(matrix[1,:])-2): #column
###            currentage = matrix[i,j]
###            nextage = matrix[i,j+1]
###
###            if currentage > nextage:
###                print('next age wrong!')
###                print('i',i)
###                print('j',j)
###                
###
###test_convert_to_single_particle(age_subset)
###    
###    
###    
###
###with open(parent_dir + '\data\output\simulations\Simulation_' + savename + '.dictionary', 'wb') as config_dictionary_file:
###        
###         pickle.dump(trajectory_data, config_dictionary_file)
##
##
##
#









