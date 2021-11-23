# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 13:50:01 2021

@author: quint
"""

import numpy as np



Coastgrids = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData\Coastgrids.npy')

gridnumber = 0
GridNumberMask = np.zeros((len(Coastgrids[:,0]), len(Coastgrids[0,:])))


for row in range(len(Coastgrids[:,0])):
    for column in range(len(Coastgrids[0,:])):
        
        grid = Coastgrids[row,column]
        
        if grid == 1:
            
            gridnumber = gridnumber + 1
            
            GridNumberMask[row,column] = gridnumber
        
np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Data\Output\Galapagos_RGEMS_FieldData/GridNumberMask.npy', GridNumberMask)

