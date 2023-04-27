import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df_max = pd.read_csv('C:/Users/91956/OneDrive - students.iiserpune.ac.in/sims/glu_200_2/max_values/max_values.csv')
Cacyto = df_max.Max_Ca_cytoplasm.tolist()
Cacyto.sort()
print(Cacyto)

bins = list(np.arange(0.5,20,0.5))
plt.xlabel('Calcium peak (micromolar)')
plt.ylabel('Number of Calcium events')
plt.title('Peak amplitude histogram (continuous vs train stimulus)')
plt.hist(Cacyto, bins=bins, edgecolor = 'black')