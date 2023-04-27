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
path = os.getcwd()
csv_files = glob.glob(os.path.join('C:/Users/91956/OneDrive - students.iiserpune.ac.in/Sem_06/Semester Project_Prof. Suhita/AD_PMCA_1/*.csv'))

#list containing max values of correspoding concentrations from all the files in a given simulation run

Cacyto = []
CaER = []
IP3conc = []
x = 0
threshold = 0.3 #for calicum events to be considered (Anup's paper)

# loop over the list of csv files
for f in csv_files:
    # read the csv file
    df = pd.read_csv(f)
    cacyto = df['Ca_cyto'] 
    Cacyto_elements = [x for x in cacyto if x > threshold]
    Cacyto.append(Cacyto_elements)
    #Cacyto.append(max(df['Ca_cyto']))
    #CaER.append(max(df['Ca_ER']))
    #IP3conc.append(max(df['IP3_conc']))
    x = x + 1
    print(x/50)

#dictionary of lists
dict = {'Max_Ca_cytoplasm': Cacyto} #'Max_Ca_ER': CaER, 'Max_IP3': IP3conc
df2 = pd.DataFrame(dict)

df2.to_csv('C:/Users/91956/OneDrive - students.iiserpune.ac.in/Sem_06/Semester Project_Prof. Suhita/AD_PMCA_1/ca_event/peak_amplitudes.csv')

'''
Cacyto.sort() #sort the lists in ascending order
CaER.sort()
IP3conc.sort()

bins = list(np.arange(0.5,17,0.5))
plt.hist(Cacyto, bins=bins, edgecolor = 'black')
'''
