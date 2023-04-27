import pandas as pd
import os
import glob
import numpy as np
from scipy.signal import find_peaks, peak_widths
import matplotlib.pyplot as plt

rise_times = []

for filename in os.listdir('C:/Users/91956/OneDrive - students.iiserpune.ac.in/sims/glu_200_2'): 
        if filename.endswith(".csv"):
            file_directory = os.path.join('C:/Users/91956/OneDrive - students.iiserpune.ac.in/sims/glu_200_2', filename)
            df = pd.read_csv(file_directory)
            Cacytosignal = df.Ca_cyto
            peak_indices, _ = find_peaks(Cacytosignal, height = 0.3) #idk how but adding the _ after peak variable gives only peak indices, 
            #otherwise the output contains both the peak indices and peak heights
            results_half = peak_widths(Cacytosignal, peak_indices, rel_height=0.5)
            results_20 = peak_widths(Cacytosignal, peak_indices, rel_height=0.2)
            results_80 = peak_widths(Cacytosignal, peak_indices, rel_height=0.8)
            
            for i, peak in enumerate(peak_indices):
                rise_time = results_20[2][i] - results_80[2][i]
                rise_times.append(1000*rise_time)

dict = {'rise_times': rise_times}
df2 = pd.DataFrame(dict)
df2.to_csv('C:/Users/91956/OneDrive - students.iiserpune.ac.in/sims/glu_200_2/max_values/rise_times.csv')


'''
bins = list(np.arange(0,210,10))
plt.hist(rise_times, bins=bins, edgecolor = 'black')
'''
'''
peak_width_80pc = peak_widths(Cacytosignal, peaks, rel_height=0.2)
peak_width_20pc = peak_widths(Cacytosignal, peaks, rel_height=0.8)
print(peak_width_20pc)
'''