# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 10:10:57 2023

@author: 91956
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import peak_widths, find_peaks
import os
import glob


#good ol' raw plot for concentrations
'''
z = []
q = []
for i in range(0,301,10):
    z.append(i)
for j in range(0,401):
    q.append(j)    
glu400_30 = pd.read_csv('C:\semester project\AD_PMCA_1/output0.csv')
plt.plot(glu400_30.Time, glu400_30.Ca_cyto, label = 'Ca_cytoplasm')
#plt.plot(glu400_30.Time, glu400_30.Ca_ER, label = 'Ca_ER')
plt.plot(glu400_30.Time, glu400_30.IP3_conc, label = 'IP3_conc')
plt.xlabel('Time (in seconds)')
plt.ylabel('Concentration (in ${\mu}$M)')
plt.title('Calcium transients (AD PMCA with single cluster containing 10000 IP$_{3}$R channels)')
#plt.xticks(z)
#plt.yticks(q)
#plt.yscale('log')
plt.legend(loc='upper right')
#plt.figure(figsize=(50,50), dpi=500)
plt.show()
'''


j = -1
for filename in os.listdir('C:/Users/91956/OneDrive - students.iiserpune.ac.in/Sem_06/Semester Project_Prof. Suhita/control_1'): 
        if filename.endswith(".csv"):
            file_directory = os.path.join('C:/Users/91956/OneDrive - students.iiserpune.ac.in/Sem_06/Semester Project_Prof. Suhita/control_1', filename)
            AD_PMCA_i = pd.read_csv(file_directory)
            ca_cyto = AD_PMCA_i.Ca_cyto
            #print (ca_cyto)
            peak_indices, _ = find_peaks(ca_cyto, height = 0.3, distance = 100) #indices for calcium events, the distance of 100 corresponds to 5ms
            time = peak_indices*0.00005 #multiplying to convert it into seconds                     
            peak_values = ca_cyto[peak_indices]
            dict = {'time': time, 'peak_values': peak_values}
            df = pd.DataFrame(dict)
            j = j + 1
            df.to_csv('C:/Users/91956/OneDrive - students.iiserpune.ac.in/Sem_06/Semester Project_Prof. Suhita/control_1/peaks/peak' + str(j) + '.csv')
          


'''
plt.scatter(df.time, df.peak_values, label = 'Ca_cyto')
plt.xlabel('Time (in seconds)')
plt.ylabel('Concentration (in ${\mu}$M)')
plt.xticks(np.arange(0, 105, 5))
plt.legend(loc = 'upper right')
plt.show()
'''
