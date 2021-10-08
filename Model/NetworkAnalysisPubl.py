# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 10:42:20 2021

@author: quint
"""


import sys
sys.path.append(r"C:\Users\quint\Documents\Quinten_studie\Publicatie\Modules")
import numpy as np
import pandas as pd
from NetworkAnalysisModule import NetworkAnalysis as Nw

bordercurrentlist = [0, 0.0000003, 0.0000006,0.0000009,0.0000012]

scaling = 5000

savename = 'bc is,' + str(bordercurrentlist[0])
TransitionMatrix = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\TransitionMatrices\TM' + savename + '.npy')
GridsDataFrame = pd.read_csv(r'C:\Users\quint\Documents\Quinten_studie\Master_Thesis\Data\Output_Data\ReleaseLocations\GridsDataFrame.csv')

NAOBC = Nw(TransitionMatrix, GridsDataFrame)
NAOBC.CreateDirectedGraph()
NAOBC.BetwneennesCentrality()
NAOBC.ConstructNetworkDataFrame()
NAOBC.PlotCentrality(scaling = scaling,  title = str(bordercurrentlist[0]))



savename = 'bc is,' + str(bordercurrentlist[1])
TransitionMatrix = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\TransitionMatrices\TM' + savename + '.npy')
GridsDataFrame = pd.read_csv(r'C:\Users\quint\Documents\Quinten_studie\Master_Thesis\Data\Output_Data\ReleaseLocations\GridsDataFrame.csv')

NAO = Nw(TransitionMatrix, GridsDataFrame)
NAO.CreateDirectedGraph()
NAO.BetwneennesCentrality()
NAO.ConstructNetworkDataFrame()
NAO.PlotCentrality(scaling = scaling, title = str(bordercurrentlist[1]))

savename = 'bc is,' + str(bordercurrentlist[2])
TransitionMatrix = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\TransitionMatrices\TM' + savename + '.npy')
GridsDataFrame = pd.read_csv(r'C:\Users\quint\Documents\Quinten_studie\Master_Thesis\Data\Output_Data\ReleaseLocations\GridsDataFrame.csv')

NAO = Nw(TransitionMatrix, GridsDataFrame)
NAO.CreateDirectedGraph()
NAO.BetwneennesCentrality()
NAO.ConstructNetworkDataFrame()
NAO.PlotCentrality(scaling = scaling, title = str(bordercurrentlist[2]))

savename = 'bc is,' + str(bordercurrentlist[3])
TransitionMatrix = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\TransitionMatrices\TM' + savename + '.npy')
GridsDataFrame = pd.read_csv(r'C:\Users\quint\Documents\Quinten_studie\Master_Thesis\Data\Output_Data\ReleaseLocations\GridsDataFrame.csv')

NAO = Nw(TransitionMatrix, GridsDataFrame)
NAO.CreateDirectedGraph()
NAO.BetwneennesCentrality()
NAO.ConstructNetworkDataFrame()
NAO.PlotCentrality(scaling = scaling, title = str(bordercurrentlist[3]))

savename = 'bc is,' + str(bordercurrentlist[4])
TransitionMatrix = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\TransitionMatrices\TM' + savename + '.npy')
GridsDataFrame = pd.read_csv(r'C:\Users\quint\Documents\Quinten_studie\Master_Thesis\Data\Output_Data\ReleaseLocations\GridsDataFrame.csv')

NAO = Nw(TransitionMatrix, GridsDataFrame)
NAO.CreateDirectedGraph()
NAO.BetwneennesCentrality()
NAO.ConstructNetworkDataFrame()
NAO.PlotCentrality(scaling = scaling, title = str(bordercurrentlist[4]))