import os
# import vtk
import numpy as np
from Anatomic_Features import Anatomy_Features
# import pydicom
import pandas as pd
from scipy import interpolate
from dicompylercore import dicomparser,dvh, dvhcalc
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.svm import SVR
import matplotlib.pyplot as plt


defalt_path = "C:\\Users\\xhuae08006\\OneDrive - Elekta\\Desktop\\demo\\NPCDATA_DVH_DTH\\"
strt_name = 'Eye R'
pt_id = ['002','004','005','006','009','0010','0015','0016','0017','0019','0020','0023','0024','0025','0026','0028']
### DATA Collection ###
# DVH
Vol = pd.read_csv(os.path.join(defalt_path,"Volume.csv"))
X = Anatomy_Features(1)
Vol.index = pt_id
test_path = 'E:\\GitFolder\\RL-Application-in-TPS\\AUTO-PLANNING\\AutoTemplateTuning\\projects\\DVH_prediction\\DATA_'
DATA = X.DVH_data_extration(test_path,pt_id,structure_name=strt_name)
DATA.to_csv(os.path.join(defalt_path,"DVH_"+strt_name+"_17.csv"))

# DTH
# Target_Name = 'PCTV'
# OAR_Name = 'parotid l
DTH_STAT_DATA = []
for item in pt_id:
    rtssfile = "E:/download20200908/XHTOMODATA RAW/"+item+"/"+item+"_StrctrSets.dcm"
    rtss = dicomparser.DicomParser(rtssfile)
    rtstructures = rtss.GetStructures()
    for i in rtstructures.keys():
        if 'pctv' in rtstructures[i]['name'].lower():
            Target_ID = rtstructures[i]['id']
        elif strt_name.lower() in rtstructures[i]['name'].lower():
            OAR_ID = rtstructures[i]['id']
    print(item)
    DTH_data = X._calculate_DTH(rtssfile,Target_ID, OAR_ID, wplt=False)
    DTH_STAT_DATA.append(DTH_data['DTH'])

y_new = X.DTH_data_extration(DTH_STAT_DATA,pt_id)
y_new.to_csv(os.path.join(defalt_path,"DTH_"+ strt_name +"_17.csv"))

## Just try to calculate DTH all the cases
DTH,DVH = {},{}
log_inf = []
import os
import numpy as np
from Anatomic_Features import Anatomy_Features
from dicompylercore import dicomparser,dvh, dvhcalc
X = Anatomy_Features(1)
path = "D:\\NBIA_HNSCC_DATA\\HNSCC_rtss\\"
dose_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_rtdose\\"
log_path = "C:\\Users\\xhuae08006\\OneDrive - Elekta\\Desktop\\demo\\"
filenames = os.listdir(path)
for name in filenames:
    ## calculate Parotid DTH
    DTH[name.split('.')[0]] = {}
    rtss = dicomparser.DicomParser(os.path.join(path,name))
    rtstructures = rtss.GetStructures()
    Target_ID,OAR_ID = {},{}
    for i in rtstructures.keys():
        if 'ptv' in rtstructures[i]['name'].lower():
            Target_ID[rtstructures[i]['name']] = rtstructures[i]['id']
        elif 'parotid' in rtstructures[i]['name'].lower() and 'sub' not in rtstructures[i]['name'].lower():
            OAR_ID[rtstructures[i]['name']] = rtstructures[i]['id']
    # for key1 in Target_ID.keys():
    #     for key2 in OAR_ID.keys():
    #         print("PT {} for DTH is calculated between {} and {}".format(name.split('.')[0],key1,key2))
    #         DTH_data = X._calculate_DTH(os.path.join(path,name),Target_ID[key1], OAR_ID[key2], wplt=False)
    #         DTH[name.split('.')[0]][key1+"_VS_"+key2] = DTH_data
    #         if DTH_data['DTH'] != []:
    #             if not np.all(DTH_data['DTH'][0] > 0):
    #                 log_inf.append("PT"+name.split('.')[0]+key1+"_VS_"+key2+"DTH\n")
    ## calculate Parotid DVH
    DVH[name.split('.')[0]] = {}
    for key3 in OAR_ID.keys():
        print("PT {} is calculating {}'s DVH".format(name.split('_')[0]+"rtdose",key3))
        try:
            calcdvh = dvhcalc.get_dvh(os.path.join(path,name), os.path.join(dose_path,name.split('_')[0]+"_rtdose.dcm"), OAR_ID[key3])
            DVH[name.split('.')[0]][key3] = {'Volume':calcdvh.counts/calcdvh.counts[0]*100, 'Dose':calcdvh.bins }
        except AttributeError:
            print("'Dataset' object has no attribute 'ContourGeometricType'")
            DVH[name.split('.')[0]][key3] = {}

DVH_Left_Parotid = {}
DVH_Right_Parotid = {}

for key in DVH.keys():
    for key2 in DVH[key].keys():
        if 'parotid' in key2.lower() and 'lt' in key2.lower():
            DVH_Left_Parotid[key] = DVH[key][key2]
        elif 'parotid' in key2.lower() and 'l' in key2.lower():
            DVH_Left_Parotid[key] = DVH[key][key2]
        elif 'parotid' in key2.lower() and 'rt' in key2.lower():
            DVH_Right_Parotid[key] = DVH[key][key2]
        elif 'parotid' in key2.lower() and 'r ' in key2.lower() or ' r' in key2.lower() or 'right' in key2.lower():
            DVH_Right_Parotid[key] = DVH[key][key2]
        else:
            print(key2)


import pandas as pd
DTH = pd.DataFrame(DTH)
DTH.to_csv("C:\\Users\\xhuae08006\\OneDrive - Elekta\\Desktop\\demo\\NBIA_Head_Neck_DATA\\DTH.csv")
DVH = pd.DataFrame(DVH)
DVH.to_csv("C:\\Users\\xhuae08006\\OneDrive - Elekta\\Desktop\\demo\\NBIA_Head_Neck_DATA\\DVH.csv")




## Create a clean rtdose dicom for further DVH calculation ##
input_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx\\"
output_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_rtdose\\"
clean_rtss_name = os.listdir(input_path)
for name in clean_rtss_name:
    dcm = pydicom.read_file(os.path.join("D:\\NBIA_HNSCC_DATA\\HNSCC_rtdose\\",name+'_rtdose.dcm'))
    dcm.save_as(os.path.join(output_path,name+'_rtdose.dcm'))


## calculate DVH ##
## Just try to calculate DTH all the cases
DVH = {}
log_inf = []
import os
import numpy as np
from Anatomic_Features import Anatomy_Features
from dicompylercore import dicomparser,dvh, dvhcalc
# X = Anatomy_Features(1)
rtss_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_rtss\\"
dose_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_rtdose\\"
log_path = "C:\\Users\\xhuae08006\\OneDrive - Elekta\\Desktop\\demo\\"
filenames = os.listdir(rtss_path)
for name in filenames:
    ## calculate Parotid DTH
    # DTH[name.split('.')[0]] = {}
    rtss = dicomparser.DicomParser(os.path.join(rtss_path,name))
    rtstructures = rtss.GetStructures()
    Target_ID,OAR_ID = {},{}
    for i in rtstructures.keys():
        if 'ptv' in rtstructures[i]['name'].lower():
            Target_ID[rtstructures[i]['name']] = rtstructures[i]['id']
        elif 'BrainStem' in rtstructures[i]['name'] or "SpinalCord" in rtstructures[i]['name'] and 'sub' not in rtstructures[i]['name'].lower():
            OAR_ID[rtstructures[i]['name']] = rtstructures[i]['id']
    # for key1 in Target_ID.keys():
    #     for key2 in OAR_ID.keys():
    #         print("PT {} for DTH is calculated between {} and {}".format(name.split('.')[0],key1,key2))
    #         DTH_data = X._calculate_DTH(os.path.join(path,name),Target_ID[key1], OAR_ID[key2], wplt=False)
    #         DTH[name.split('.')[0]][key1+"_VS_"+key2] = DTH_data
    #         if DTH_data['DTH'] != []:
    #             if not np.all(DTH_data['DTH'][0] > 0):
    #                 log_inf.append("PT"+name.split('.')[0]+key1+"_VS_"+key2+"DTH\n")

    ## calculate Parotid DVH
    DVH[name.split('.')[0]] = {}
    for key3 in OAR_ID.keys():
        print("PT {} is calculating {}'s DVH".format(name.split('_')[0]+"rtdose",key3))
        try:
            calcdvh = dvhcalc.get_dvh(os.path.join(rtss_path,name), os.path.join(dose_path,name.split('_')[0]+"_rtdose.dcm"), OAR_ID[key3])
            DVH[name.split('.')[0]][key3] = {'Volume':calcdvh.counts/calcdvh.counts[0]*100, 'Dose':calcdvh.bins }
        except AttributeError:
            print("'Dataset' object has no attribute 'ContourGeometricType'")
            DVH[name.split('.')[0]][key3] = {}

DVH_Left_Parotid = {}
DVH_Right_Parotid = {}

for key in DVH.keys():
    for key2 in DVH[key].keys():
        if  key2 == 'BrainStem':
            DVH_Left_Parotid[key] = DVH[key][key2]
        elif key2 == 'SpinalCord':
            DVH_Right_Parotid[key] = DVH[key][key2]


Left_Parotid_Dose_max =[(DVH_Left_Parotid[key]['Dose'][-1],key) for key in DVH_Left_Parotid.keys()]
Right_Parotid_Dose_max = [(DVH_Right_Parotid[key]['Dose'][-1],key) for key in DVH_Right_Parotid.keys()]
Left_Parotid_Dose_Index = np.arange(0,max(Left_Parotid_Dose_max)[0],0.01)
Right_Parotid_Dose_Index = np.arange(0,max(Right_Parotid_Dose_max)[0],0.01)

Left_Parotid_DVH_Matrix = np.zeros([len(DVH_Left_Parotid.keys()),Left_Parotid_Dose_Index.shape[0]]) 
Right_Parotid_DVH_Matrix = np.zeros([len(DVH_Right_Parotid.keys()),Right_Parotid_Dose_Index.shape[0]]) 
columns_left,columns_right = [],[]
for i,item in enumerate(DVH_Left_Parotid.keys()):
    Left_Parotid_DVH_Matrix[i,0:DVH_Left_Parotid[item]['Volume'].shape[0]] = DVH_Left_Parotid[item]['Volume']
    columns_left.append(item)
for i,item in enumerate(DVH_Right_Parotid.keys()):
    Right_Parotid_DVH_Matrix[i,0:DVH_Right_Parotid[item]['Volume'].shape[0]] = DVH_Right_Parotid[item]['Volume']
    columns_right.append(item)

import pandas as pd
Left_Parotid_DVH_Matrix_ = pd.DataFrame(Left_Parotid_DVH_Matrix.T,columns=columns_left,index = Left_Parotid_Dose_Index)
Right_Parotid_DVH_Matrix_ = pd.DataFrame(Right_Parotid_DVH_Matrix.T,columns=columns_right,index = Right_Parotid_Dose_Index)
Left_Parotid_DVH_Matrix_.to_csv("C:\\Users\\xhuae08006\\OneDrive - Elekta\\Desktop\\demo\\NBIA_Head_Neck_DATA\\DVH_BrainStem_Clean.csv")
Right_Parotid_DVH_Matrix_.to_csv("C:\\Users\\xhuae08006\\OneDrive - Elekta\\Desktop\\demo\\NBIA_Head_Neck_DATA\\DVH_SpinalCord_Clean.csv")


import matplotlib.pyplot as plt
fig = plt.figure
Right_Parotid_DVH_Matrix_.plot()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('Dose(cGy)')
plt.ylabel('Volume(%)')
plt.title('DVH_Parotid R')
plt.grid(b=True)
# axes.set_ylim([0,100])
plt.show()


import matplotlib.pyplot as plt
fig = plt.figure
Left_Parotid_DVH_Matrix_.plot()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('Dose(cGy)')
plt.ylabel('Volume(%)')
plt.title('DVH_Parotid L')
plt.grid(b=True)
# axes.set_ylim([0,100])
plt.show()


## To Calculate DTH of each parotid with PTV ##
## Just try to calculate DTH all the cases
DTH_Left,DTH_Right = {},{}
left_log_inf,right_log_inf = [],[]
import os
import numpy as np
from Anatomic_Features import Anatomy_Features
from dicompylercore import dicomparser,dvh, dvhcalc
X = Anatomy_Features(1)
rtss_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_rtss\\"
dose_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_rtdose\\"
log_path = "C:\\Users\\xhuae08006\\OneDrive - Elekta\\Desktop\\demo\\"
filenames = os.listdir(rtss_path)
for name in filenames:
    ## calculate Parotid DTH
    DTH_Left[name.split('.')[0]] = {}
    DTH_Right[name.split('.')[0]] = {}
    rtss = dicomparser.DicomParser(os.path.join(rtss_path,name))
    rtstructures = rtss.GetStructures()
    Target_ID,OAR_ID = {},{}
    for i in rtstructures.keys():
        if 'ptv' in rtstructures[i]['name'].lower():
            Target_ID[rtstructures[i]['name']] = rtstructures[i]['id']
        elif 'BrainStem' in rtstructures[i]['name'] or "SpinalCord" in rtstructures[i]['name'] and 'sub' not in rtstructures[i]['name'].lower():
            OAR_ID[rtstructures[i]['name']] = rtstructures[i]['id']
    for key1 in Target_ID.keys():
        for key2 in OAR_ID.keys():
            if key2 == "SpinalCord":
                print("PT {} for DTH is calculated between {} and {}".format(name.split('.')[0],key1,key2))
                DTH_data = X._calculate_DTH(os.path.join(rtss_path,name),Target_ID[key1], OAR_ID[key2], wplt=False)
                DTH_Left[name.split('.')[0]][key1+"_VS_"+key2] = DTH_data
                if DTH_data['DTH'] != []:
                    if not np.all(DTH_data['DTH'][0] >= 0):
                        left_log_inf.append("PT"+name.split('.')[0]+key1+"_VS_"+key2+"DTH\n")
            # elif key2 == "Parotid_R":
            #     print("PT {} for DTH is calculated between {} and {}".format(name.split('.')[0],key1,key2))
            #     DTH_data = X._calculate_DTH(os.path.join(rtss_path,name),Target_ID[key1], OAR_ID[key2], wplt=False)
            #     DTH_Right[name.split('.')[0]][key1+"_VS_"+key2] = DTH_data
            #     if DTH_data['DTH'] != []:
            #         if not np.all(DTH_data['DTH'][0] >= 0):
            #             right_log_inf.append("PT"+name.split('.')[0]+key1+"_VS_"+key2+"DTH\n")
min_index,max_index = [],[]
count = 0
for key in DTH_Left.keys():
    for key2 in DTH_Left[key].keys():
        print(key)
        if DTH_Left[key][key2]['DTH'] != []:
            min_index.append((min(DTH_Left[key][key2]['DTH'][0]),key))
            max_index.append((max(DTH_Left[key][key2]['DTH'][0]),key))
            count += 1

Left_Parotid_Index = np.arange(min(min_index)[0],max(max_index)[0],0.01)
Left_Matrix = np.zeros([count,Left_Parotid_Index.shape[0]])
from scipy import interp
count1 = 0
columns = []
for key in DTH_Left.keys():
    for key2 in DTH_Left[key].keys():
        print(key)
        if DTH_Left[key][key2]['DTH'] != []:
            y_new = interp(Left_Parotid_Index,DTH_Left[key][key2]['DTH'][0],DTH_Left[key][key2]['DTH'][1])
            Left_Matrix[count1,0:y_new.shape[0]] = y_new
            count1 += 1
            columns.append(key+'__'+key2)

import pandas as pd
DTH_Left_Matrix = pd.DataFrame(Left_Matrix.T,index = Left_Parotid_Index,columns = columns)
DTH_Left_Matrix.to_csv("C:\\Users\\xhuae08006\\OneDrive - Elekta\\Desktop\\demo\\NBIA_Head_Neck_DATA\\DTH_SpinalCord_Clean.csv")



from Anatomic_Features import Standard_Tools
data_folder = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_PTVs\\"
Y = Standard_Tools(data_folder)

Y.Rename_STRTNAME_in_rtssdcm()

#DTH_new_columns = []
#
#for name_DVH in ML_Model.DVH_Data.columns:
#    
#    os.mkdir(os.path.join(DVH_DTH_data_path,name_DVH+"_Fig"))    
#    #DVH show
#    print('Figure!')
#    fig = plt.figure()
#    # axes= fig.add_axes([0.1,0.1,0.8,0.8])
#    ML_Model.DVH_Data[name_DVH].plot()
#    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
#    plt.xlabel('Dose(Gy)')
#    plt.ylabel('volume(%)')
#    plt.title('DVH_ParotidL')
#    # axes.set_ylim([0,1])
#    # plt.xlim(-11.61,35.34)
#    plt.grid(b=True)
#    fig_name1 = os.path.join(DVH_DTH_data_path,name_DVH+"_Fig",name_DVH+"_DVH.tif")
#    plt.savefig(fig_name1,dpi=200, bbox_inches='tight')
#    plt.show()
#    
#    name_DTH = [item for item in ML_Model.DTH_Data.columns if name_DVH in item]
#    print('Figure!')
#    fig = plt.figure()
#    # axes= fig.add_axes([0.1,0.1,0.8,0.8])
#    ML_Model.DTH_Data[name_DTH].plot()
#    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
#    plt.xlabel('distance(mm)')
#    plt.ylabel('volume(%)')
#    plt.title('DTH_ParotidL')
#    # axes.set_ylim([0,1])
#    # plt.xlim(-11.61,35.34)
#    plt.grid(b=True)
#    fig_name = os.path.join(DVH_DTH_data_path,name_DVH+"_Fig",name_DVH+"_DTH.tif")
#    plt.savefig(fig_name,dpi=200, bbox_inches='tight')
#    plt.show()
#
#    ## further extract DTH
#    SS = np.abs(ML_Model.DTH_Data[name_DTH].index - 0)
#    indx = np.where(SS == np.min(SS))
#    Temp = ML_Model.DTH_Data[name_DTH].iloc[indx[0][0],:]
#    indx_1 = np.where(Temp == np.max(Temp))
#    DTH_new_columns.append(Temp.index[indx_1[0][0]])


###################### Rename SpinalCord and BrainStem ########################
# import pydicom
# import os
# input_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_rtss\\"
# # output_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_rtss\\"
# dirty_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_rtss\\"
# clean_rtss_name = os.listdir(input_path)
# for name in clean_rtss_name:
#     print("Name:{}".format(name))
#     dcm = pydicom.read_file(os.path.join(input_path,name))
#     count = 0
#     for ROI in dcm.StructureSetROISequence:
#         if 'expanded' not in ROI.ROIName.lower() and \
#                'expd' not in ROI.ROIName.lower() and \
#                 'exp' not in ROI.ROIName.lower() and \
#                 'hot' not in ROI.ROIName.lower() and \
#              'expand' not in ROI.ROIName.lower() and \
#                 'new' not in ROI.ROIName.lower() and \
#                 'prv' not in ROI.ROIName.lower() and \
#               '0.5cm' not in ROI.ROIName.lower() and \
#               'cords' not in ROI.ROIName.lower() and \
#                  'no' not in ROI.ROIName.lower() and \
#                  'mp' not in ROI.ROIName.lower() and \
#                 'pcv' not in ROI.ROIName.lower() and \
#               'avoid' not in ROI.ROIName.lower() and \
#               '0.5mm' not in ROI.ROIName.lower() and \
#            'planning' not in ROI.ROIName.lower() and \
#            ' 2' not in ROI.ROIName.lower():
#             if 'cord' in ROI.ROIName.lower():
#                 print('original name:{}'.format(ROI.ROIName)) 
#                 dcm.StructureSetROISequence[dcm.StructureSetROISequence.index(ROI)].ROIName = 'SpinalCord'
#                 count +=1
#             if 'stem' in ROI.ROIName.lower():
#                 print('original name:{}'.format(ROI.ROIName))
#                 dcm.StructureSetROISequence[dcm.StructureSetROISequence.index(ROI)].ROIName = 'BrainStem'
#                 count +=1
#     print("Number:{}".format(count))
#     dcm.save_as(os.path.join(input_path,name))

# ######################### Standardization of HNSCC parotid left and right names ####################################
# import pydicom
# import os
# input_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx\\"
# output_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_rtss\\"
# dirty_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_rtss\\"
# clean_rtss_name = os.listdir(input_path)
# for name in clean_rtss_name:
#     print("Name:{}".format(name))
#     dcm = pydicom.read_file(os.path.join(dirty_path,name+'_rtss.dcm'))
#     for ROI in dcm.StructureSetROISequence:
#         if 'parotid' in ROI.ROIName.lower() and 'sub' not in ROI.ROIName.lower() and \
#          'avoid' not in ROI.ROIName.lower() and 'low' not in ROI.ROIName.lower() and \
#            '20' not in ROI.ROIName.lower() and 'tail' not in ROI.ROIName.lower() and \
#            'sud' not in ROI.ROIName.lower() and 'ptv' not in ROI.ROIName.lower() and \
#            'parotid2' not in ROI.ROIName.lower() and 'parotids' not in ROI.ROIName.lower() and \
#              'push' not in ROI.ROIName.lower() and 'hot' not in ROI.ROIName:
#             if 'lt' in ROI.ROIName.lower() or 'left' in ROI.ROIName.lower() or 'L' in ROI.ROIName or 'l ' in ROI.ROIName or 'lft' in ROI.ROIName:
#                 print('original name:{}'.format(ROI.ROIName))
#                 dcm.StructureSetROISequence[dcm.StructureSetROISequence.index(ROI)].ROIName = 'Parotid_L'
#                 print(dcm.StructureSetROISequence[dcm.StructureSetROISequence.index(ROI)].ROIName)
#             elif 'rt' in ROI.ROIName.lower() or 'right' in ROI.ROIName.lower() or 'R' in ROI.ROIName or 'r ' in ROI.ROIName:
#                 print('original name:{}'.format(ROI.ROIName))
#                 dcm.StructureSetROISequence[dcm.StructureSetROISequence.index(ROI)].ROIName = 'Parotid_R'
#                 print(dcm.StructureSetROISequence[dcm.StructureSetROISequence.index(ROI)].ROIName)
    
#     dcm.save_as(os.path.join(output_path,name+'_rtss.dcm'))            