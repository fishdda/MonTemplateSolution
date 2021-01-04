'''
   1. To display cost function in DVH to vividly illustrate the cost functions' effective range in DVH curves
   2. To help physicist describe the inner principle of cost functions
'''

## extract DVH data from a sample DICOM
## calculate serial value in each dose points
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
DVH_path = "C:\\Users\\xhuae08006\\OneDrive - Elekta\\Desktop\\demo\\NBIA_Head_Neck_DATA\\DVH_BrainStem_Clean.csv"
DVH_Stem = pd.read_csv(DVH_path)
Sample_DVH = DVH_Stem['HNSCC-01-0003_rtss']

D_0 = 20   # unit:Gy
K   = 10   # Power Law
X = np.arange(0,Sample_DVH.shape[0]/100,0.01)
f = pd.DataFrame((X/D_0)**K)
Sample_DVH_Plot = pd.concat([Sample_DVH,f],axis=1)
Sample_DVH_Plot.columns = ['HNSCC-01-0002_rtss','K=10&D_0=20Gy']
plt.figure()
Sample_DVH_Plot.plot()

plt.xlabel("Dose(Gy)")
plt.ylabel("Volume(%)")
plt.grid()
plt.show()



# Parallel function
D_0 = 30   # unit:Gy
K   = 15   # Power Law
X = np.arange(0,Sample_DVH.shape[0]/100,0.01)
MM = 1/(1+(D_0/X)**K)
MM_1 = MM/np.max(MM)*100
KK = 5*(MM_1[1:] - MM_1[0:-1])/0.01

f = pd.DataFrame(KK)
Sample_DVH_Plot = pd.concat([Sample_DVH_Plot,f],axis=1)
Sample_DVH_Plot.columns = ['HNSCC-01-0002_rtss','K=5&D_0=25Gy','K=10&D_0=25Gy','K=1&D_0=25Gy','K=5&D_0=30Gy','K=15&D_0=30Gy']
plt.figure()
Sample_DVH_Plot.plot()

plt.xlabel("Dose(Gy)")
plt.ylabel("Volume(%)")
plt.grid()
plt.savefig("C:\\Users\\xhuae08006\\OneDrive - Elekta\\Desktop\\CFs_display.png",dpi=200, bbox_inches='tight')
plt.show()