# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 18:00:34 2020

@author: xhuae08006
"""

## to run Anatomic_Features' PCA_SVR modules ##

from Anatomic_Features import PCA_SVR_Model
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd

DVH_DTH_data_path = "E:\\demo\\NBIA_Head_Neck_DATA\\" 
strt_name         = "BrainStem"
#epsilon           = 0.06
epsilon_1         = np.arange(0,1,0.01)
epsilon_2         = np.arange(0,1,0.01)
epsilon_3         = np.arange(0,1,0.01)
epsilon_4         = np.arange(0,1,0.01)
epsilon_5         = np.arange(0,1,0.01)
ML_Model = PCA_SVR_Model(DVH_DTH_data_path,strt_name)

ML_Model.PCA_Compression(num_components=5,wplt=False,DVH_file_name="DVH_"+strt_name+"_Clean.csv",DTH_file_name="DTH_"+strt_name+"_Clean.csv")

ML_Model.Training_VS_Validation_Set(rtss_Volume_path=os.path.join(DVH_DTH_data_path,"Volume_Clean_New.csv"),structure_name=strt_name)

#y_test_pd, y_pred_pd = ML_Model.SVR_Correlation_Multi_output_version(epsilon)

y_test_pd, y_pred_pd = ML_Model.SVR_Correlation_Single_output_version(epsilon_1,epsilon_2,epsilon_3,epsilon_4,epsilon_5)
## to see the difference between PCA compressed data and uncompressed data ###
#BackData = ML_Model.SC_DVH.inverse_transform(np.dot(ML_Model.principalComponents_DVH,ML_Model.PCA_DVH.components_))




DVH_Valid_Pred, DVH_Valid_Truth = ML_Model.DVH_Reconstruct()


Final_Data = pd.concat([DVH_Valid_Pred,DVH_Valid_Truth],axis=1,sort=True)

#Final_Data.plot()
name = list(set([item.split('_PCA_SVR')[0] for item in Final_Data.columns]))
Test1 = Final_Data[[name[0]+'_PCA_SVR_Pred',name[0]]]
Test2 = Final_Data[[name[1]+'_PCA_SVR_Pred',name[1]]]
Test3 = Final_Data[[name[2]+'_PCA_SVR_Pred',name[2]]]
Test1['difference'] = np.abs(Test1.iloc[:,0]- Test1.iloc[:,1])
Test2['difference'] = np.abs(Test2.iloc[:,0]- Test2.iloc[:,1])
Test3['difference'] = np.abs(Test3.iloc[:,0]- Test3.iloc[:,1])

fig = plt.figure()                
Test1.plot()
plt.xlabel('Dose(Gy)')
plt.ylabel('Volume(%)')
plt.title('DVH_'+strt_name)
plt.grid(b=True)
plt.savefig(name[0] +"_DVH.tif",dpi=200, bbox_inches='tight')
plt.show()

fig = plt.figure()                
Test2.plot()
plt.xlabel('Dose(Gy)')
plt.ylabel('Volume(%)')
plt.title('DVH_'+strt_name)
plt.grid(b=True)
plt.savefig(name[1] + "_DVH.tif",dpi=200, bbox_inches='tight')
plt.show()

fig = plt.figure()                
Test3.plot()
plt.xlabel('Dose(Gy)')
plt.ylabel('Volume(%)')
plt.title('DVH_'+strt_name)
plt.grid(b=True)
plt.savefig(name[2] + "_DVH.tif",dpi=200, bbox_inches='tight')
plt.show()



## plot the data distribution
SS = ML_Model.DVH_Data.T
DATA = SS.describe()
DATA = DATA.T
DATA_PLOT = DATA[['mean','min','max','25%','50%','75%']]
fig = plt.figure()                
DATA_PLOT.plot()
plt.xlabel('Dose(Gy)')
plt.ylabel('Volume(%)')
plt.title('DVH_Statistics_'+strt_name)
plt.grid(b=True)
plt.savefig(strt_name + "_DVHDataStatistics.tif",dpi=200, bbox_inches='tight')
plt.show()

