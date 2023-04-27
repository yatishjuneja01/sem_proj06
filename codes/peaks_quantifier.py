# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 22:55:34 2023

@author: 91956
"""

import pandas as pd
import os
import glob
import numpy as np
from scipy.signal import find_peaks, peak_widths
import matplotlib.pyplot as plt
peak_amplitudes = []
times = []
j = -1
for filename in os.listdir('C:/Users/91956/OneDrive - students.iiserpune.ac.in/Sem_06/Semester Project_Prof. Suhita/control_200/peaks'): 
        if filename.endswith(".csv"):
            file_directory = os.path.join('C:/Users/91956/OneDrive - students.iiserpune.ac.in/Sem_06/Semester Project_Prof. Suhita/control_200/peaks', filename)
            AD_PMCA_i = pd.read_csv(file_directory)
            PEAK = AD_PMCA_i.peak_values
            TIME = AD_PMCA_i.time
            peak_amplitudes.append(PEAK)
            times.append(TIME)
            #print (ca_cyto)
            #peak_values = ca_cyto[peak_indices]
            #dict = {'time': time, 'peak_values': peak_values}
            #df = pd.DataFrame(dict)
            #j = j + 1
            #df.to_csv('C:/Users/91956/OneDrive - students.iiserpune.ac.in/Sem_06/Semester Project_Prof. Suhita/AD_PMCA_1/peaks/peak' + str(j) + '.csv')
            
# Concatenate the arrays into a single array
peak_amplitudes = np.concatenate(peak_amplitudes).ravel()
times = np.concatenate(times).ravel()

'''
bins = list(np.arange(0,100,1))
plt.ylabel('Frequency')
plt.xlabel('Time (in seconds)')
plt.title('Calcium event occurence frequency histogram (control with 10000 IP${_3}$Rs receptors)')
plt.hist(times, bins=bins, edgecolor = 'black')
plt.xticks(np.arange(0, 105, 5))
plt.show()
'''

bins = list(np.arange(0.3,21,0.3))
plt.ylabel('Frequency')
plt.xlabel('Peak amplitude (in ${\mu}$M)')
plt.title('Peak Amplitude Histogram (control with 200 IP${_3}$R receptors)')
plt.hist(peak_amplitudes, bins=bins, edgecolor = 'black')
plt.xticks(np.arange(0.3, 21, 0.9))
plt.show()


'''
bins = list(np.arange(0,10,1))
plt.ylabel('Frequency')
plt.xlabel('Peak amplitude (in  ${\mu}$M)')
plt.title('Peak Amplitude Histogram (AD PMCA with 1 IP${_3}$R receptor)')
plt.hist(peak_amplitudes, bins=bins, edgecolor = 'black')
plt.show()
'''