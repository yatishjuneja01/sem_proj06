# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 15:20:35 2023

@author: 91956
"""

import pandas as pd
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

#list containing max values of correspoding concentrations from all the files in a given simulation run
Cacyto = []
CaER = []
IP3conc = []
x = 0

# loop over the list of csv files
for filename in os.listdir('/home/user/Downloads/sem_proj_yatish/sims/glu_200_30/'):
    if filename.endswith(".csv"):
        file_directory = os.path.join('/home/user/Downloads/sem_proj_yatish/sims/glu_200_30/', filename)
        df = pd.read_csv(file_directory)
        Cacyto.append(max(df['Ca_cyto']))
        CaER.append(max(df['Ca_ER']))
        IP3conc.append(max(df['IP3_conc']))
        x = x + 1
        print(x/2)

#dictionary of lists
dict = {'Max_Ca_cytoplasm': Cacyto, 'Max_Ca_ER': CaER, 'Max_IP3': IP3conc}
df2 = pd.DataFrame(dict)

df2.to_csv('/home/user/Downloads/sem_proj_yatish/sims/glu_200_30/max_values/max_values.csv')


Cacyto.sort() #sort the lists in ascending order
CaER.sort()
IP3conc.sort()

'''
bins = list(np.arange(0.5,20,0.5))
plt.xlabel('Calcium peak (micromolar)')
plt.ylabel('Number of Calcium events')
plt.title('Peak amplitude histogram for train stimulus (200)')
plt.hist(Cacyto, bins=bins, edgecolor = 'black')
'''
