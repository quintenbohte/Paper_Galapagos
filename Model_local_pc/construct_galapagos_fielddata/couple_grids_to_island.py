# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 12:17:33 2021

@author: quint
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

coastgrids = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Model_local_pc\data\input\galapagos_field_data\Coastgrids.npy')
gridsdataframe = pd.read_csv(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Model_local_pc\data\input\galapagos_field_data\GridsDataFrame.csv')
gridnumber = np.load(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Model_local_pc\data\input\galapagos_field_data\GridNumberMask.npy')


island1 = [361,360,359,356,353,352,357,349,348,345,344,341,340,339,
           336,330,331,324,325,326,327,332,333,334,328,322,315,308,300,301,
           294,285,293,284,292,291,283,282,275,274,273,262,250,237,230,220,221,
           203,204,205,206,207,208,209,210,222,231,232,238,251,264,263,276,286,
           252,253,239,240,233,223,211,212,195,196,188,189,190,180,181,176,175,161,
           160,174,173,172,159,158,157,156,145,144,134,133,118,117,102,101,91,90,81,
           76,68,55,56,50,51,40,41,42,43,44,45,46,35,36,37,38,39,47,48,49,52,53,54,
           57,58,59,60,69,70,71,77,82,92,103,104,119,135,147,146,163,162,177,185,184,183,
           182,191,192,197,213,224,234,243,242,241,254,267,266,265,278,277,288,287,
           295,304,303,302,309,316,323,329,335,338,337,343,342,347,346,350,355,354,358, 272,351]

island1.sort()
island2 = [8,9,10,11,12,16,17,24,26,29,30,34,33,32,31,28,27,25,23,15]

island3 = [1,2,3,4,7,14,22,21,20,19,18,13,5,6]

island4= [249,248,247,235,227,226,225,218,217,216,215,214,199,198,200,201,219,228,229,
          193,186,178,164,148,136,137,138,120,121,122,105,106,107,
          108,109,123,124,125,139,140,150,165,179,187,194,202,236,149]

island5 = [83,84,85,95,111,129,128,127,126,110,93,94]

island6 = [171,170,169,154,168,167,166,153,152,151,142,141,130,
           114,113,97,96,88,87,86,78,72,73,61,62,63,64,65,66,67,
           74,75,79,80,89,98,99,100,112,115,116,131,132,143,155]

island7 = [321,320,319,318,317,310,305,297,296,289,279,280,268,
           269,270,255,256,257,258,259,244,245,246,260,261,271,281,290,
           299,298,307,306,314,313,312,311]

island8 = [380,381,382,383,385,387,390,389,392,393,391,388,386,384]


island9 = [379,378,377,376,375,372,370,367,368,362,363,364, 365, 366, 369, 371,374,373]

allgridsisland = island1+island2+island3+island4+island5+island6+island7+island8+island9

allgridsunique = set(allgridsisland)

total_island_list = []
total_island_list.append(island1)
total_island_list.append(island2)
total_island_list.append(island3)
total_island_list.append(island4)
total_island_list.append(island5)
total_island_list.append(island6)
total_island_list.append(island7)
total_island_list.append(island8)
total_island_list.append(island9)


island_number = []



contains = 68 in island1
print(contains)

for grid in gridsdataframe['Grid_Number']:
    for island in range(len(total_island_list)):
        islandlist = total_island_list[island]
        print(islandlist)
        contain = grid in islandlist
        if contain:
            island_number.append(island+1)
            break
        
    

gridsdataframe['island_number'] = island_number

###########################MAKE ISLAND NUMBER MASK###################################

df=gridsdataframe

island_number = np.asarray(gridsdataframe['island_number'])


island_number_mask = np.zeros((  len(gridnumber[:,1]),len(gridnumber[1,:])))

for i in range(len(gridnumber[:,1])): #vert
    for j in range(len(gridnumber[1,:])):
        
        grid = int(gridnumber[i,j])
        
        if grid != 0:
        
            island_number_grid = island_number[grid-1]
            island_number_mask[i,j] = island_number_grid
            


gridsdataframe.to_csv(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Model_local_pc\data\input\galapagos_field_data\GridsDataFrame.csv', index = False)
np.save(r'C:\Users\quint\Documents\Quinten_studie\Publicatie\Model_local_pc\data\input\galapagos_field_data\islandnumber_mask', island_number_mask)





















