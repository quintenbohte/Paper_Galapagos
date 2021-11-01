# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 09:13:53 2021

@author: quint
"""

import numpy as np
from parcels import Field, FieldSet, JITParticle, ScipyParticle 
from parcels import ParticleFile, ParticleSet, Variable, VectorField, ErrorCode
from parcels.tools.converters import GeographicPolar 
from datetime import timedelta as delta
from os import path
from glob import glob
import numpy as np
import dask
import math
import xarray as xr
from netCDF4 import Dataset
import warnings
import matplotlib.pyplot as plt
import pickle
warnings.simplefilter('ignore', category=xr.SerializationWarning)
from operator import attrgetter
import pickle


class Simulation:
    
    def __init__(self, length_simulation, advection_duration, 
                 output_frequency, repeatdt, 
                 deltatime, savename, 
                 domain, path_to_velocity, 
                 path_to_grid, data_in, bordercurrent = None):
        
        self.ReleaseLat = list(np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\ReleaseLocations\ReleaseLat.npy'))
        self.ReleaseLon = list(np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\ReleaseLocations\ReleaseLon.npy'))
        self.length_simulation = length_simulation
        self.advection_duration = advection_duration
        self.output_frequency = output_frequency
        self.domain = domain
        self.savename = savename
        self.path_to_velocity = path_to_velocity
        self.path_to_grid = path_to_grid
        self.data_in = data_in
        self.repeatdt = repeatdt
        self.deltatime = deltatime
        self.bordercurrent = bordercurrent
        
        self.fieldset = None
        self.londomain = None
        self.latdomain = None
        self.loncor = None
        self.latcor = None
        self.lon = None
        self.lat = None
        self.iy_min = None
        self.ix_min = None
        self.iy_max = None
        self.ix_max = None
        self.pset = None
        self.outfile = None
        self.kernels = None
        self.trajectory_data = None
        self.beached_particles = None
        
        
    def getclosest_ij(self, lats,lons,latpt,lonpt):    
        
        dist_lat = (lats-latpt)**2          # find squared distance of every point on grid
        dist_lon = (lons-lonpt)**2
        minindex_lat = dist_lat.argmin()    # 1D index of minimum dist_sq element
        minindex_lon = dist_lon.argmin()
        
        return minindex_lat, minindex_lon   # Get 2D index for latvals and lonvals arrays from 1D index
        
    def create_fieldset(self):
        
        dfile = Dataset(self.data_in+'\RGEMS3_Surf_grid.nc')
        self.lon = dfile.variables['XC'][:]
        self.lat = dfile.variables['YC'][:]
        self.loncor = dfile.variables['XG'][:]
        self.latcor = dfile.variables['YG'][:]
        self.iy_min, self.ix_min = self.getclosest_ij(self.lat, self.lon, self.domain[2], self.domain[0])
        self.iy_max, self.ix_max = self.getclosest_ij(self.lat, self.lon, self.domain[3], self.domain[1])
        
        self.londomain = self.lon[self.ix_min:self.ix_max]
        self.latdomain = self.lat[self.iy_min:self.iy_max]
        
        ############ CREATE FIELDSET ####################
        
        PathToVelocity = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input\RGEMS_2010_Surf.nc'
        PathToGrid = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input\RGEMS3_Surf_grid.nc'
        
        VelocityField = sorted(glob(PathToVelocity))
        GridField = glob(PathToGrid)
        
        files_MITgcm = {'U': {'lon': GridField, 'lat': GridField, 'data': VelocityField},                                             
                    'V': {'lon': GridField, 'lat': GridField, 'data': VelocityField}}     
        
        variables = {'U': 'UVEL', 'V': 'VVEL'}
        dimensions = {'lon': 'XG', 'lat': 'YG', 'time': 'time'}
        
        indices = {'lon': range(self.ix_min,self.ix_max), 
                   'lat': range(self.iy_min,self.iy_max)}
        
        
        self.fieldset = FieldSet.from_mitgcm(files_MITgcm,
                                             variables,
                                             dimensions,
                                             indices = indices)
        

    def add_fields(self):
        
        fieldset = self.fieldset
        
        landmask = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\Landmask.npy')
        coastgrids = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\Coastgrids.npy')
        coastgrids = np.roll(coastgrids,1,0)
        gridnumbermask = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\GridNumberMask.npy')
        gridnumbermask = np.roll(gridnumbermask,1,0)
        
        fieldset.add_constant('advection_duration',self.advection_duration)    
        fieldset.add_constant('lon_max',self.lon[self.ix_max])
        fieldset.add_constant('lon_min',self.lon[self.ix_min])
        fieldset.add_constant('lat_max',self.lat[self.iy_max])
        fieldset.add_constant('lat_min',self.lat[self.iy_min])
        fieldset.add_constant('bordercurrent',self.bordercurrent)


        
        fieldset.add_field(Field('landmask',
                                 data = landmask,
                                 lon = self.londomain,
                                 lat = self.latdomain,
                                 mesh='spherical',
                                 interp_method = 'nearest'))
        
        fieldset.add_field(Field('coastgrids',
                                 data = coastgrids,
                                 lon = self.londomain,
                                 lat = self.latdomain,
                                 mesh='spherical',
                                 interp_method = 'nearest'))
        
        fieldset.add_field(Field('gridnumbermask',
                                 data = gridnumbermask,
                                 lon = self.londomain,
                                 lat = self.latdomain,
                                 mesh='spherical',
                                 interp_method = 'nearest'))
        
        
        
        BorU = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\borderCurrentU.npy')
        BorUlon = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\BorderUlon.npy')
        BorUlat = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\BorderUlat.npy')
        
        BorV = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\borderCurrentV.npy')
        BorVlon = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\BorderVlon.npy')
        BorVlat = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\BorderVlat.npy')
        
        fieldset.add_field(Field('BorderCurrentU', BorU, BorUlon, BorUlat))
        fieldset.add_field(Field('BorderCurrentV', BorV, BorVlon, BorVlat))
            
        
        DistanceFromShore = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\distance2shore.npy')
        x = np.linspace(self.domain[0],self.domain[1],len(self.londomain))
        y = np.linspace(self.domain[2],self.domain[3],len(self.londomain))
        lon, lat = np.meshgrid(x,y)
        fieldset.add_field(Field('distance2shore', DistanceFromShore, lon, lat))
        



                
    def AdvectionRK4(particle, fieldset, time):
        """ Only advect particles that are not out of bounds"""
        if (particle.lon < fieldset.lon_max and
            particle.lon > fieldset.lon_min and
            particle.lat < fieldset.lat_max and
            particle.lat > fieldset.lat_min):
        
            (u1, v1) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
            lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
            
            if (lon1 < fieldset.lon_max and
                lon1 > fieldset.lon_min and
                lat1 < fieldset.lat_max and
                lat1 > fieldset.lat_min):
    
                (u2, v2) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat1, lon1]
                lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)
    
                if (lon2 < fieldset.lon_max and
                    lon2 > fieldset.lon_min and
                    lat2 < fieldset.lat_max and
                    lat2 > fieldset.lat_min):
    
                    (u3, v3) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat2, lon2]
                    lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)
    
                    if (lon3 < fieldset.lon_max and
                        lon3 > fieldset.lon_min and
                        lat3 < fieldset.lat_max and
                        lat3 > fieldset.lat_min):
    
                        (u4, v4) = fieldset.UV[time + particle.dt, particle.depth, lat3, lon3]
                        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
                        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt


    def create_particle(self):
        
        class GalapagosParticle(JITParticle):
            distance = Variable('distance', initial = 0) 
            age = Variable('age', dtype=np.float32, initial = 0.)
            beached = Variable('beached', dtype=np.int32, initial = 0.) 
            coastcell = Variable('coastcell', dtype=np.int32, initial = 0.)
            gridnumber = Variable('gridnumber', dtype=np.int32, initial = 0.)
            distancetoshore = Variable('distancetoshore', dtype=np.int32, initial = 0.)
#    StrengthBorderCurrent = Variable('StrengthBorderCurrent', dtype=np.int32, initial = StrengthBorderCurrent)
    
    
    
        if self.repeatdt == None:
            print('NoRepeat')
            self.pset = ParticleSet(fieldset = self.fieldset,
                                   pclass = GalapagosParticle,
                                   lon = self.ReleaseLon,
                                   lat = self.ReleaseLat,
                                   time = 0)
        else:
            print('repeat')
            self.pset = ParticleSet(fieldset = self.fieldset,
                                   pclass = GalapagosParticle,
                                   lon = self.ReleaseLon,
                                   lat = self.ReleaseLat,
                                   repeatdt = delta(hours = self.repeatdt),
                                   time = 0)  
            
    
    def add_kernels(self):
    
        
                     
        def AdvectionRK4(particle, fieldset, time):
            """ Only advect particles that are not out of bounds"""
            if (particle.lon < fieldset.lon_max and
                particle.lon > fieldset.lon_min and
                particle.lat < fieldset.lat_max and
                particle.lat > fieldset.lat_min):
            
                (u1, v1) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
                lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
                
                if (lon1 < fieldset.lon_max and
                    lon1 > fieldset.lon_min and
                    lat1 < fieldset.lat_max and
                    lat1 > fieldset.lat_min):
        
                    (u2, v2) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat1, lon1]
                    lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)
        
                    if (lon2 < fieldset.lon_max and
                        lon2 > fieldset.lon_min and
                        lat2 < fieldset.lat_max and
                        lat2 > fieldset.lat_min):
        
                        (u3, v3) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat2, lon2]
                        lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)
        
                        if (lon3 < fieldset.lon_max and
                            lon3 > fieldset.lon_min and
                            lat3 < fieldset.lat_max and
                            lat3 > fieldset.lat_min):
        
                            (u4, v4) = fieldset.UV[time + particle.dt, particle.depth, lat3, lon3]
                            particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
                            particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
            
            
        def Age(fieldset, particle, time):
            """ Delete particles when reaching age specified by advection_duration """
            particle.age = particle.age + math.fabs(particle.dt)
            if particle.age > fieldset.advection_duration*86400:
                particle.delete()
        
        def SampleInfo(fieldset, particle, time):
            if (particle.lon < fieldset.lon_max and
                particle.lon > fieldset.lon_min and
                particle.lat < fieldset.lat_max and
                particle.lat > fieldset.lat_min):
                
                particle.distance = fieldset.distance[time, particle.depth, particle.lat, particle.lon]
                particle.island = fieldset.island[time, particle.depth, particle.lat, particle.lon]    
        
        def DeleteParticle(particle, fieldset, time):
            if particle.lon < fieldset.lon_min+0.1 or particle.lon > fieldset.lon_max-0.1 or particle.lat < fieldset.lat_min+0.1 or particle.lat > fieldset.lat_max-0.1: 
                particle.delete()
            
        
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
                
        def AntiBeachNudging(particle,fieldset,time):
            """    
            The nudging current is 1 m s^-1, which ought to be sufficient to overpower
            any coastal current (I hope) and push our particle back out to sea so as to
            not get stuck
            
            update 11/03/2020: Following tests and discussions with Cleo, the nudging 
            current will now kick in starting at 500m from the coast, since otherwise 
            the particles tended to get stuck if we used the velocity treshhold. 
            """
            particle.distancetoshore = fieldset.distance2shore[time,particle.depth,particle.lat,particle.lon]
            if particle.distancetoshore < 2:
                print('test')
                borUab,borVab=fieldset.BorderCurrentU[time, particle.depth, particle.lat, particle.lon],fieldset.BorderCurrentV[time, particle.depth, particle.lat, particle.lon]
                particle.lon += fieldset.bordercurrent*borUab*particle.dt
                particle.lat += fieldset.bordercurrent*borVab*particle.dt
                
                
                
        self.outfile = self.pset.ParticleFile(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Simulations\Simulation_' + self.savename + '.nc', outputdt = delta(hours = self.output_frequency))
    
        self.kernels = (self.pset.Kernel(AdvectionRK4) + self.pset.Kernel(Age) +
                           self.pset.Kernel(beachtesting) + self.pset.Kernel(coasttesting) +
                           self.pset.Kernel(DeleteParticle) + self.pset.Kernel(gridnumbertesting))

    
    
    def execute(self):
        self.pset.execute(self.kernels,
                     runtime=delta(days=self.length_simulation),
                     dt=delta(hours=self.deltatime),
                     output_file=self.outfile)
        
        self.outfile.export()
        self.outfile.close()
            
    
    def add_trajectory_data(self):
        
        self.trajectory_data = xr.open_dataset(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Simulations\Simulation_' + self.savename + '.nc')
    
    def determine_beached_particles(self):
        
        beached = self.trajectory_data['beached'].data
        beached = np.nan_to_num(beached)
        AllParticles = np.sum(beached, axis = 1)
        BeachedBarticles = np.count_nonzero(AllParticles)
        self.beached_particles = BeachedBarticles/len(AllParticles)
        
        

        
    def to_dict(self):
        dictionary = {'length_simulation':self.length_simulation,
                      'advection_duration':self.advection_duration,
                      'output_frequency':self.output_frequency,
                      'domain':self.domain,
                      'savename':self.savename,
                      'repeatdt':self.repeatdt,
                      'trajectory_data':self.trajectory_data,
                      'deltatime': self.deltatime, 
                      'londomain':self.londomain, 
                      'latdomain':self.latdomain, 
                      'loncor':self.loncor,
                      'latcor':self.latcor, 
                      'iy_min':self.iy_min,
                      'ix_min': self.ix_min, 
                      'iy_max': self.iy_max, 
                      'ix_max': self.ix_max,
                      'beached_particles': self.beached_particles}
        
        

        with open(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Simulations\dict' +self.savename+  '.dictionary', 'wb') as config_dictionary_file:
        
          pickle.dump(dictionary, config_dictionary_file)
            

    
    def run_simulation(self):
        self.create_fieldset()
        self.add_fields()
        self.create_particle()
        self.add_kernels()
        self.execute()
        self.add_trajectory_data()
        self.determine_beached_particles()
        self.to_dict()
        
        





class Simulation_set_back:
    
    def __init__(self, length_simulation, advection_duration, 
                 output_frequency, repeatdt, 
                 deltatime, savename, 
                 domain, path_to_velocity, 
                 path_to_grid, data_in, bordercurrent = None):
        
        self.ReleaseLat = list(np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\ReleaseLocations\ReleaseLat.npy'))
        self.ReleaseLon = list(np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\ReleaseLocations\ReleaseLon.npy'))
        self.length_simulation = length_simulation
        self.advection_duration = advection_duration
        self.output_frequency = output_frequency
        self.domain = domain
        self.savename = savename
        self.path_to_velocity = path_to_velocity
        self.path_to_grid = path_to_grid
        self.data_in = data_in
        self.repeatdt = repeatdt
        self.deltatime = deltatime
        self.bordercurrent = bordercurrent
        
        self.fieldset = None
        self.londomain = None
        self.latdomain = None
        self.loncor = None
        self.latcor = None
        self.lon = None
        self.lat = None
        self.iy_min = None
        self.ix_min = None
        self.iy_max = None
        self.ix_max = None
        self.pset = None
        self.outfile = None
        self.kernels = None
        self.trajectory_data = None
        self.beached_particles = None
        
        
    def getclosest_ij(self, lats,lons,latpt,lonpt):    
        
        dist_lat = (lats-latpt)**2          # find squared distance of every point on grid
        dist_lon = (lons-lonpt)**2
        minindex_lat = dist_lat.argmin()    # 1D index of minimum dist_sq element
        minindex_lon = dist_lon.argmin()
        
        return minindex_lat, minindex_lon   # Get 2D index for latvals and lonvals arrays from 1D index
        
    def create_fieldset(self):
        
        dfile = Dataset(self.data_in+'\RGEMS3_Surf_grid.nc')
        self.lon = dfile.variables['XC'][:]
        self.lat = dfile.variables['YC'][:]
        self.loncor = dfile.variables['XG'][:]
        self.latcor = dfile.variables['YG'][:]
        self.iy_min, self.ix_min = self.getclosest_ij(self.lat, self.lon, self.domain[2], self.domain[0])
        self.iy_max, self.ix_max = self.getclosest_ij(self.lat, self.lon, self.domain[3], self.domain[1])
        
        self.londomain = self.lon[self.ix_min:self.ix_max]
        self.latdomain = self.lat[self.iy_min:self.iy_max]
        
        ############ CREATE FIELDSET ####################
        
        PathToVelocity = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input\RGEMS_2010_Surf.nc'
        PathToGrid = r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Input\RGEMS3_Surf_grid.nc'
        
        VelocityField = sorted(glob(PathToVelocity))
        GridField = glob(PathToGrid)
        
        files_MITgcm = {'U': {'lon': GridField, 'lat': GridField, 'data': VelocityField},                                             
                    'V': {'lon': GridField, 'lat': GridField, 'data': VelocityField}}     
        
        variables = {'U': 'UVEL', 'V': 'VVEL'}
        dimensions = {'lon': 'XG', 'lat': 'YG', 'time': 'time'}
        
        indices = {'lon': range(self.ix_min,self.ix_max), 
                   'lat': range(self.iy_min,self.iy_max)}
        
        
        self.fieldset = FieldSet.from_mitgcm(files_MITgcm,
                                             variables,
                                             dimensions,
                                             indices = indices)
        

    def add_fields(self):
        
        fieldset = self.fieldset
        
        landmask = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\Landmask.npy')
        coastgrids = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\Coastgrids.npy')
        coastgrids = np.roll(coastgrids,1,0)
        gridnumbermask = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\GridNumberMask.npy')
        gridnumbermask = np.roll(gridnumbermask,1,0)
        
        fieldset.add_constant('advection_duration',self.advection_duration)    
        fieldset.add_constant('lon_max',self.lon[self.ix_max])
        fieldset.add_constant('lon_min',self.lon[self.ix_min])
        fieldset.add_constant('lat_max',self.lat[self.iy_max])
        fieldset.add_constant('lat_min',self.lat[self.iy_min])
        fieldset.add_constant('bordercurrent',self.bordercurrent)


        
        fieldset.add_field(Field('landmask',
                                 data = landmask,
                                 lon = self.londomain,
                                 lat = self.latdomain,
                                 mesh='spherical',
                                 interp_method = 'nearest'))
        
        fieldset.add_field(Field('coastgrids',
                                 data = coastgrids,
                                 lon = self.londomain,
                                 lat = self.latdomain,
                                 mesh='spherical',
                                 interp_method = 'nearest'))
        
        fieldset.add_field(Field('gridnumbermask',
                                 data = gridnumbermask,
                                 lon = self.londomain,
                                 lat = self.latdomain,
                                 mesh='spherical',
                                 interp_method = 'nearest'))
        
        
        
        BorU = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\borderCurrentU.npy')
        BorUlon = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\BorderUlon.npy')
        BorUlat = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\BorderUlat.npy')
        
        BorV = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\borderCurrentV.npy')
        BorVlon = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\BorderVlon.npy')
        BorVlat = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\BorderVlat.npy')
        
        fieldset.add_field(Field('BorderCurrentU', BorU, BorUlon, BorUlat))
        fieldset.add_field(Field('BorderCurrentV', BorV, BorVlon, BorVlat))
            
        
        DistanceFromShore = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\distance2shore.npy')
        x = np.linspace(self.domain[0],self.domain[1],len(self.londomain))
        y = np.linspace(self.domain[2],self.domain[3],len(self.londomain))
        lon, lat = np.meshgrid(x,y)
        fieldset.add_field(Field('distance2shore', DistanceFromShore, lon, lat))
        



                
    def AdvectionRK4(particle, fieldset, time):
        """ Only advect particles that are not out of bounds"""
        if (particle.lon < fieldset.lon_max and
            particle.lon > fieldset.lon_min and
            particle.lat < fieldset.lat_max and
            particle.lat > fieldset.lat_min):
        
            (u1, v1) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
            lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
            
            if (lon1 < fieldset.lon_max and
                lon1 > fieldset.lon_min and
                lat1 < fieldset.lat_max and
                lat1 > fieldset.lat_min):
    
                (u2, v2) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat1, lon1]
                lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)
    
                if (lon2 < fieldset.lon_max and
                    lon2 > fieldset.lon_min and
                    lat2 < fieldset.lat_max and
                    lat2 > fieldset.lat_min):
    
                    (u3, v3) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat2, lon2]
                    lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)
    
                    if (lon3 < fieldset.lon_max and
                        lon3 > fieldset.lon_min and
                        lat3 < fieldset.lat_max and
                        lat3 > fieldset.lat_min):
    
                        (u4, v4) = fieldset.UV[time + particle.dt, particle.depth, lat3, lon3]
                        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
                        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt





    def create_particle(self):
        
        class GalapagosParticle(JITParticle):
            distance = Variable('distance', initial = 0) 
            age = Variable('age', dtype=np.float32, initial = 0.)
            beached = Variable('beached', dtype=np.int32, initial = 0.) 
            coastcell = Variable('coastcell', dtype=np.int32, initial = 0.)
            gridnumber = Variable('gridnumber', dtype=np.int32, initial = 0.)
            distancetoshore = Variable('distancetoshore', dtype=np.int32, initial = 0.)
            startlon = Variable('startlon', dtype=np.float32, initial = self.ReleaseLon)
            startlat = Variable('startlat', dtype=np.float32, initial = self.ReleaseLat)
            delta_time = Variable('delta_time', dtype=np.float32, initial = self.deltatime)



        if self.repeatdt == None:
            print('NoRepeat')
            self.pset = ParticleSet(fieldset = self.fieldset,
                                   pclass = GalapagosParticle,
                                   lon = self.ReleaseLon,
                                   lat = self.ReleaseLat,
                                   time = 0)
        else:
            print('repeat')
            self.pset = ParticleSet(fieldset = self.fieldset,
                                   pclass = GalapagosParticle,
                                   lon = self.ReleaseLon,
                                   lat = self.ReleaseLat,
                                   repeatdt = delta(hours = self.repeatdt),
                                   time = 0)  
            
    
    def add_kernels(self):
    
        
                     
        def AdvectionRK4(particle, fieldset, time):
            """ Only advect particles that are not out of bounds"""
            if (particle.lon < fieldset.lon_max and
                particle.lon > fieldset.lon_min and
                particle.lat < fieldset.lat_max and
                particle.lat > fieldset.lat_min):
            
                (u1, v1) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
                lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
                
                if (lon1 < fieldset.lon_max and
                    lon1 > fieldset.lon_min and
                    lat1 < fieldset.lat_max and
                    lat1 > fieldset.lat_min):
        
                    (u2, v2) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat1, lon1]
                    lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)
        
                    if (lon2 < fieldset.lon_max and
                        lon2 > fieldset.lon_min and
                        lat2 < fieldset.lat_max and
                        lat2 > fieldset.lat_min):
        
                        (u3, v3) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat2, lon2]
                        lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)
        
                        if (lon3 < fieldset.lon_max and
                            lon3 > fieldset.lon_min and
                            lat3 < fieldset.lat_max and
                            lat3 > fieldset.lat_min):
        
                            (u4, v4) = fieldset.UV[time + particle.dt, particle.depth, lat3, lon3]
                            particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
                            particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
            
            
        def Age(fieldset, particle, time):
            """ Delete particles when reaching age specified by advection_duration """
            particle.age = particle.age + math.fabs(particle.dt)
#            if particle.age > fieldset.advection_duration*86400:
#                particle.delete()
        
        def SampleInfo(fieldset, particle, time):
            if (particle.lon < fieldset.lon_max and
                particle.lon > fieldset.lon_min and
                particle.lat < fieldset.lat_max and
                particle.lat > fieldset.lat_min):
                
                particle.distance = fieldset.distance[time, particle.depth, particle.lat, particle.lon]
                particle.island = fieldset.island[time, particle.depth, particle.lat, particle.lon]    
        
#        def DeleteParticle(particle, fieldset, time):
#            if particle.lon < fieldset.lon_min+0.1 or particle.lon > fieldset.lon_max-0.1 or particle.lat < fieldset.lat_min+0.1 or particle.lat > fieldset.lat_max-0.1: 
#                particle.delete()
#            
#        
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
                
        def AntiBeachNudging(particle,fieldset,time):
            """    
            The nudging current is 1 m s^-1, which ought to be sufficient to overpower
            any coastal current (I hope) and push our particle back out to sea so as to
            not get stuck
            
            update 11/03/2020: Following tests and discussions with Cleo, the nudging 
            current will now kick in starting at 500m from the coast, since otherwise 
            the particles tended to get stuck if we used the velocity treshhold. 
            """
            particle.distancetoshore = fieldset.distance2shore[time,particle.depth,particle.lat,particle.lon]
            if particle.distancetoshore < 2:
                print('test')
                borUab,borVab=fieldset.BorderCurrentU[time, particle.depth, particle.lat, particle.lon],fieldset.BorderCurrentV[time, particle.depth, particle.lat, particle.lon]
                particle.lon += fieldset.bordercurrent*borUab*particle.dt
                particle.lat += fieldset.bordercurrent*borVab*particle.dt
                
                
        def SetParticleBack(particle,fieldset,time):
    
            if particle.age > fieldset.advection_duration*86400:
                
#                print('test')
                particle.lon = particle.startlon
                particle.lat = particle.startlat
                particle.age = 0
        
        
                
                
        self.outfile = self.pset.ParticleFile(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Simulations\Simulation_' + self.savename + '.nc', outputdt = delta(hours = self.output_frequency))
    
        self.kernels = (self.pset.Kernel(AdvectionRK4) + self.pset.Kernel(Age) +
                           self.pset.Kernel(beachtesting) + self.pset.Kernel(coasttesting) +
                            self.pset.Kernel(gridnumbertesting) + self.pset.kernel(SetParticleBack))

    
    
    def execute(self):
        self.pset.execute(self.kernels,
                     runtime=delta(days=self.advection_duration),
                     dt=delta(hours=self.deltatime),
                     output_file=self.outfile)
        
        self.outfile.export()
        self.outfile.close()
        
        self.pset.repeatdt = None
        
        
        self.pset.execute(self.kernels,
             runtime=delta(days=self.length_simulation - self.advection_duration),
             dt=delta(hours=self.deltatime),
             output_file=self.outfile)
        
        

    
    def add_trajectory_data(self):
        
        self.trajectory_data = xr.open_dataset(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Simulations\Simulation_' + self.savename + '.nc')
    
    def determine_beached_particles(self):
        
        beached = self.trajectory_data['beached'].data
        beached = np.nan_to_num(beached)
        AllParticles = np.sum(beached, axis = 1)
        BeachedBarticles = np.count_nonzero(AllParticles)
        self.beached_particles = BeachedBarticles/len(AllParticles)
        
        

        
    def to_dict(self):
        dictionary = {'length_simulation':self.length_simulation,
                      'advection_duration':self.advection_duration,
                      'output_frequency':self.output_frequency,
                      'domain':self.domain,
                      'savename':self.savename,
                      'repeatdt':self.repeatdt,
                      'trajectory_data':self.trajectory_data,
                      'deltatime': self.deltatime, 
                      'londomain':self.londomain, 
                      'latdomain':self.latdomain, 
                      'loncor':self.loncor,
                      'latcor':self.latcor, 
                      'iy_min':self.iy_min,
                      'ix_min': self.ix_min, 
                      'iy_max': self.iy_max, 
                      'ix_max': self.ix_max,
                      'beached_particles': self.beached_particles}
        
        

        with open(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Simulations\dict' +self.savename+  '.dictionary', 'wb') as config_dictionary_file:
        
          pickle.dump(dictionary, config_dictionary_file)
            
    
    def convert_to_single_part(self):
        
        data = self.trajectory_data
        
        
        output_frequency = self.output_frequency     #hours
#        repeated_days = self.advection_duration     #the amount of days that we re-release particles
        advection_duration_hours = int(self.advection_duration*(24/output_frequency))+1 #the hours a particle id advected before it is set back
        release_locations = len(self.ReleaseLon)            #Number of particles we release every time
        advected_timesteps = len(data[1,:])             #total hours that we advected particles
        number_particles = math.ceil(advected_timesteps/advection_duration_hours)       #number of particles in one trajectory
        output_trajectory = np.zeros(((len(data[:,1])*number_particles) ,advection_duration_hours)) # the output trajectory, we put the set back particles in.
        
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
                
#        
#    return output_trajectory
#    
    
    def run_simulation(self):
        self.create_fieldset()
        self.add_fields()
        self.create_particle()
        self.add_kernels()
        self.execute()
        self.add_trajectory_data()
        self.determine_beached_particles()
        self.to_dict()
        
        









#          
#    def convert_trajectories_to_single_particles(self):
#        
#        
#        self.add_trajectory_data()
#        
#        data = self.trajectory_data
#        
#        lon = data['lon'].data
#        lat = data['lat'].data
#        coastcells = data['coastcell']        
#        beached = data['beached'].data
#        grid_number = data['gridnumber'].data
#        
#        
#        
#        
#        for trajectory_number in range(len(data)):
#    
#            trajectory = data[trajectory_number,:]
#        
#            remaining_time_steps = len(trajectory)%advection_duration
#        
#            full_particle_trajectories = trajectory[0:len(trajectory)-remaining_time_steps]
#            
#            seperate_particles = np.split(full_particle_trajectories,number_particles-1)
#            last_particle = trajectory[len(full_particle_trajectories):(len(full_particle_trajectories)+remaining_time_steps)]
#        
#            particle_set = int(math.floor(trajectory_number/release_locations))
#        
#            for particle_number in range((number_particles)): 
#                
#                if particle_number < number_particles-1:
#                
#                    print('part_set:',particle_set)
#                    indices_current_particle_set = particle_set * (release_locations*number_particles)
#                    print('current_ind:',indices_current_particle_set)
#                    index = indices_current_particle_set + ((trajectory_number%release_locations) + (release_locations*particle_number))
#                    print('index:',index)
#                    
#                    output_trajectory[index,:] = seperate_particles[particle_number]
#                
#                else:
#                    print('part_set:',particle_set)
#                    indices_current_particle_set = particle_set * (release_locations*number_particles)
#                    print('current_ind:',indices_current_particle_set)
#                    index = indices_current_particle_set + ((trajectory_number%release_locations) + (release_locations*particle_number))
#                    print('index:',index)
#                    
#                    output_trajectory[index,0:217] = last_particle
#                  
                  

#
##Sim_data = SimulationObject.trajectory_data
##disttoshore = Sim_data['distancetoshore'].data
#land = Sim_data['beached'].data
#land = np.nan_to_num(land)
#AllParticles = np.sum(land, axis = 1)
#BeachedBarticles = np.count_nonzero(AllParticles)
#BeachedPercentage = BeachedBarticles/len(AllParticles)
#print(BeachedPercentag

    