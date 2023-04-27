#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 22:49:26 2023

@author: user
"""
import pandas as pd
import os
import glob
import numpy as np
from scipy.signal import find_peaks, peak_widths
import matplotlib.pyplot as plt
df = pd.read_csv('C:/Users/91956/OneDrive - students.iiserpune.ac.in/sims/glu_200_2/max_values/rise_times.csv')
rise_times = df.rise_times.tolist()

bins = list(np.arange(0,3000,50))
plt.ylabel('Number of Calcium events')
plt.xlabel('Risetime (in ms)')
plt.title('Distribution of Calcium event rise times (continuous_100trials)')
plt.hist(rise_times, bins=bins, edgecolor = 'black')
plt.show()
