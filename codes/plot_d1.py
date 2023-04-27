import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df_max = pd.read_csv('G:\My Drive\Sem_06\Semester Project_Prof. Suhita\sims\glu_200_2\max_values\max_values.csv')
Cacyto = df_max.Max_Ca_cytoplasm.tolist()
Cacyto.sort()
print(Cacyto)

z = list(np.arange(0,30,5))
bins = list(np.arange(0.5,30,0.5))
plt.yticks(z)
plt.hist(Cacyto, bins=bins, edgecolor = 'black')