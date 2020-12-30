## deal with complexity metrics
CI_path = "C:\\Users\\xhuae08006\\OneDrive - Elekta\\Desktop\\demo\\HNSCC_complexity_metrics\\"
import os
import pandas as pd

CI = {'MCS':{},'BM':{}}
filename = os.listdir(CI_path)
for name in filename:
    table = pd.read_csv(os.path.join(CI_path,name))
    CI['MCS'][name.split(".")[0]] = (table['Modulation Complexity Score']*table["Beam MU"]/table["Beam MU"].sum()).sum()
    CI['BM'][name.split(".")[0]] = (table['Beam Modulation']*table["Beam MU"]/table["Beam MU"].sum()).sum()
    CI['PI'][name.split(".")[0]] = (table['Aperture Irregularity']*table["Beam MU"]/table["Beam MU"].sum()).sum()
    CI['OOC'][name.split(".")[0]] = (table['Out-of-circle Fractional Area']*table["Beam MU"]/table["Beam MU"].sum()).sum()
    CI['SAS2mm'][name.split(".")[0]] = (table['Small Aperture Score (2mm)']*table["Beam MU"]/table["Beam MU"].sum()).sum()
    CI['SAS5mm'][name.split(".")[0]] = (table['Small Aperture Score (5mm)']*table["Beam MU"]/table["Beam MU"].sum()).sum()
    CI['SAS10mm'][name.split(".")[0]] = (table['Small Aperture Score (10mm)']*table["Beam MU"]/table["Beam MU"].sum()).sum()
    CI['SAS20mm'][name.split(".")[0]] = (table['Small Aperture Score (20mm)']*table["Beam MU"]/table["Beam MU"].sum()).sum()

CI = pd.DataFrame(CI)
CI.to_csv(CI_path+"CI.csv")

## deal with curvatures of each PTVs
import pickle
with open('E:\\NBIA_HNSCC_DATA\\curvatures.p', 'rb') as fp:
    curvatures = pickle.load(fp)

with open('E:\\NBIA_HNSCC_DATA\\curvatures_mean.p', 'rb') as fp:
    curvatures_mean = pickle.load(fp)

with open('E:\\NBIA_HNSCC_DATA\\curvatures_std.p', 'rb') as fp:
    curvatures_std = pickle.load(fp)

curvatures_ = {key.split(".")[0]:curvatures[key] for key in curvatures.keys() if "." in key}
curvatures_mean_ = {key.split(".")[0].split("_")[0]:curvatures_mean[key] for key in curvatures_mean.keys() if "." in key}
curvatures_std_ = {key.split(".")[0].split("_")[0]:curvatures_std[key] for key in curvatures_std.keys() if "." in key}

key = list(curvatures_mean_.keys())
import pandas as pd
CI = pd.read_csv(CI_path+"CI.csv")
CI_ = CI.iloc[:,1:]
CI_.index = CI.iloc[:,0]

CI_new = CI_.drop(list(S2 - S3))

curvatures_mean_new = pd.DataFrame(curvatures_mean_,index = ['curvature_mean'])
curvatures_std_new = pd.DataFrame(curvatures_std_,index = ['curvature_std'])


curvatures_std_new = curvatures_std_new.T
curvatures_mean_new = curvatures_mean_new.T
curvatures_mean_new = curvatures_mean_new.drop(list({'HNSCC-01-0089'}))
curvatures_std_new = curvatures_std_new.drop(list({'HNSCC-01-0089'}))

data = pd.concat([CI_new,curvatures_std_new,curvatures_mean_new],axis = 1)
data.to_csv("C:\\GitFolder\\dataforcorrelation.csv")


## mockprostate
import matplotlib.pyplot as plt
import os
import numpy as np
import pickle

with open('E:/demo/geometry_exploration/mock-prostate/PTV.p', 'rb') as fp:
    mockprostate = pickle.load(fp)

with open('E:/demo/geometry_exploration/C-shape-target/Cshape.p', 'rb') as fp:
    cshape = pickle.load(fp)


fig = plt.figure() 
plt.subplot(221)               
plt.hist(mockprostate['gaussian_mesh'],bins='sqrt')
plt.title('gaussian_curvatures(mockprostate)')
plt.grid(b=True)

plt.subplot(222)
plt.hist(mockprostate['mean_mesh'],bins='sqrt')
plt.title('mean_curvatures(mockprostate)')
plt.grid(b=True)

plt.subplot(223)               
plt.hist(cshape['gaussian_curvature'],bins='sqrt')
plt.title('gaussian_curvatures(cshape)')
plt.grid(b=True)

plt.subplot(224)
plt.hist(cshape['mean_curvature'],bins='sqrt')
plt.title('mean_curvatures(cshape)')
plt.grid(b=True)

# plt.savefig(Test3.columns[0] +pic_name,dpi=200, bbox_inches='tight')
plt.show()


