# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 00:53:54 2023

@author: 91956
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import peak_widths, find_peaks
import os
import glob


j = -1
for filename in os.listdir('C:/Users/91956/OneDrive - students.iiserpune.ac.in/Sem_06/Semester Project_Prof. Suhita/control_10000'): 
        if filename.endswith(".csv"):
            file_directory = os.path.join('C:/Users/91956/OneDrive - students.iiserpune.ac.in/Sem_06/Semester Project_Prof. Suhita/control_10000', filename)
            AD_PMCA_i = pd.read_csv(file_directory)
            ca_cyto = AD_PMCA_i.Ca_cyto
            #print (ca_cyto)
            peak_indices, _ = find_peaks(ca_cyto, height = 0.3, distance = 100) #indices for calcium events, the distance of 100 corresponds to 5ms
            time = peak_indices*0.00005 #multiplying to convert it into seconds                     
            peak_values = ca_cyto[peak_indices]
            dict = {'time': time, 'peak_values': peak_values}
            df = pd.DataFrame(dict)
            j = j + 1
            df.to_csv('C:/Users/91956/OneDrive - students.iiserpune.ac.in/Sem_06/Semester Project_Prof. Suhita/control_10000/peaks/peak' + str(j) + '.csv')
            