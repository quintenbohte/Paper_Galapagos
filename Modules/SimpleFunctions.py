# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 14:43:57 2021

@author: quint
"""

import numpy as np
import pandas as pd



class Functions:
    
    pass

    @classmethod 
    def ParticleOnLand(self, OutputFile):
        land = OutputFile['beached'].data
        land = np.nan_to_num(land)
        summed = np.sum(land, axis = 1)
        
        count = 0
        for som in summed:
            if som > 0:
                count += 1
        PercentageOnLand = count/len(summed)
    
        return PercentageOnLand


    def load_dict(path_to_dict):
        import pickle
        with open(path_to_dict, 'rb') as config_dictionary_file:
         
            # Step 3
            dictionary = pickle.load(config_dictionary_file)
         
            return dictionary
    








