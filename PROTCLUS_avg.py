import matplotlib.pyplot as plt
import numpy as np
import os

os.chdir('.........')          # change to dir containing PROTCLUS analyses directories of multiple repeats

data_r1 = np.array(np.loadtxt('./PROTCLUS_r1/PROTCLUS_STEP_3/raw_data.xvg'))  # path to MD sim repeat1 (r1)
data_r2 = np.array(np.loadtxt('./PROTCLUS_r2/PROTCLUS_STEP_3/raw_data.xvg'))  #                repeat2 (r2)
data_r3 = np.array(np.loadtxt('./PROTCLUS_r3/PROTCLUS_STEP_3/raw_data.xvg'))  #                ... & so on

data_avg = (data_r1 + data_r2 + data_r3) / 3                                  # obtaining average from repeats

fig, ax = plt.subplots()
im = ax.pcolormesh(data_avg, cmap='RdPu')

plt.colorbar(im, label='Number of proteins in contact')
plt.xlabel('Frames', fontsize=13)
plt.ylabel('Protein index', fontsize=13)

fig.savefig('PROTCLUS_avg.svg', format='svg', dpi=600)
fig.savefig('PROTCLUS_avg.png', dpi=1200)

