import os
import pydicom
import numpy as np
import pandas as pd
from Anatomic_Features import Anatomy_Features
from Anatomic_Features import Standard_Tools
#################################################################################################
############################ To put all the rtss.dcm into a folder ##############################
#################################################################################################
# data_folder = "D:\\NBIA_HNSCC_DATA\\HNSCC\\"
# rtss_folder = "D:\\NBIA_HNSCC_DATA\\HNSCC_rtss\\"
# rtplan_folder = "D:\\NBIA_HNSCC_DATA\\HNSCC_rtplan\\"
# rtdose_folder = "D:\\NBIA_HNSCC_DATA\\HNSCC_rtdose\\"
# data_name = os.listdir(data_folder)
# for i,item in enumerate(data_name):
#     print("HNSCC number:{} and name: {}".format(i+1,item))
#     sub_folder = os.listdir(os.path.join(data_folder,item))
#     sub_folder_names = os.listdir(os.path.join(data_folder,item,sub_folder[0]))
#     for j,jtem in enumerate(sub_folder_names):
#         sub_sub_folder_names = os.listdir(os.path.join(data_folder,item,sub_folder[0],jtem))
#         for file_name in sub_sub_folder_names:
#             dcm = pydicom.read_file(os.path.join(data_folder,item,sub_folder[0],jtem,file_name))
#             if dcm.Modality == 'RTSTRUCT':
#                 print("DCM Modality:{}".format(dcm.Modality))
#                 dcm.save_as(os.path.join(rtss_folder,item+'_rtss.dcm'))
#             elif dcm.Modality == 'RTDOSE':
#                 print("DCM Modality:{}".format(dcm.Modality))
#                 dcm.save_as(os.path.join(rtdose_folder,item+'_rtdose.dcm'))
#             elif dcm.Modality == 'RTPLAN':
#                 print("DCM Modality:{}".format(dcm.Modality))
#                 dcm.save_as(os.path.join(rtplan_folder,item+'_rtplan.dcm'))
#             elif dcm.Modality == 'CT':
#                 continue
#X = Anatomy_Features(1)
######################################################################################################
######################### To calculate the volume of each structure ##################################
######################################################################################################
path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_rtss\\"
file_names = os.listdir(path)
Volume = {item.split('.')[0]:{} for item in file_names}
# key = "HNSCC-01-0004_rtss"
for key in Volume.keys():
    rtss = pydicom.read_file(os.path.join(path,key+".dcm"),force=True)
    structures_name = {item.ROINumber:item.ROIName for item in rtss.StructureSetROISequence}
    
    for Contour in rtss.ROIContourSequence:
        ptv = []
        try:
            if len(Contour.ContourSequence) == 1:
                vol_cc = 0
            else:    
                for item in Contour.ContourSequence:
                    try:
                        ptv.extend(item.ContourData)
                    except AttributeError:
                        print("Dataset object has no attribute 'ContourData'")
                        # print(Contour.ROIName)
                        continue
                
                ptv = np.array(ptv).reshape(int(len(ptv)/3),3)
                try:
                    vol = X._calculate_structure_volume(ptv)
                    vol_cc = vol/1000
                except IndexError:
                    print("Index 1 is out of bounds for axis 0 with size 1")
                    vol_cc = 0
        except AttributeError:
            print("Dataset' object has no attribute 'ContourSequence")
            continue
        
        Volume[key][structures_name[Contour.ReferencedROINumber]] = vol_cc
    
    print("pt:{}, and num of structures:{}".format(key,len(Volume[key].keys())))

Volume_data = pd.DataFrame(Volume)
Volume_data.to_csv(os.path.join(path,"Volume.csv"))

#################################################################################################
################# To make the data compatible with Monaco Dicom Data ############################
#################################################################################################

# data_folder = "D:\\NBIA_HNSCC_DATA\\HNSCC\\"
# destination_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_MonacoData\\"
# data_name = os.listdir(data_folder)
# for i,item in enumerate(data_name):
#     os.mkdir(os.path.join(destination_path,item))
#     print("HNSCC number:{} and name: {}".format(i+1,item))
#     sub_folder = os.listdir(os.path.join(data_folder,item))
#     sub_folder_names = os.listdir(os.path.join(data_folder,item,sub_folder[0]))
#     for j,jtem in enumerate(sub_folder_names):
#         sub_sub_folder_names = os.listdir(os.path.join(data_folder,item,sub_folder[0],jtem))
#         for file_name in sub_sub_folder_names:
#             dcm = pydicom.read_file(os.path.join(data_folder,item,sub_folder[0],jtem,file_name))
#             if dcm.Modality == 'RTSTRUCT':
#                 print("DCM Modality:{}".format(dcm.Modality))
#                 dcm.save_as(os.path.join(destination_path,item,item+'_rtss.dcm'))
#             elif dcm.Modality == 'RTDOSE':
#                 print("DCM Modality:{}".format(dcm.Modality))
#                 dcm.save_as(os.path.join(destination_path,item,item+'_rtdose.dcm'))
#             elif dcm.Modality == 'RTPLAN':
#                 print("DCM Modality:{}".format(dcm.Modality))
#                 dcm.save_as(os.path.join(destination_path,item,item+'_rtplan.dcm'))
#             elif dcm.Modality == 'CT':
#                 dcm.save_as(os.path.join(destination_path,item,file_name))


#################################################################################################
####################### To clean out the patient ID with parotid  ###############################
#################################################################################################
import pandas as pd
import os
import pydicom
rtss_dicom_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx\\"
path = "C:\\Users\\xhuae08006\\OneDrive - Elekta\\Desktop\\demo\\NBIA_Head_Neck_DATA\\"
file_names = os.listdir(rtss_dicom_path)
structures_name = {}
for i,item in enumerate(file_names):
    dicom_names = os.listdir(os.path.join(rtss_dicom_path,item))
    for dcm_name in dicom_names:
        rtss = pydicom.read_file(os.path.join(rtss_dicom_path,item,dcm_name),force=True)
        if rtss.Modality == 'RTSTRUCT':
            structures_name[item.split('.')[0]] = {item.ROINumber:item.ROIName for item in rtss.StructureSetROISequence if 'ptv' in item.ROIName.lower()}

Struct_data = pd.DataFrame(structures_name)
Struct_data.to_csv(os.path.join(path,"clean_structure_names_ptvs.csv"))

structures_name_clean = {} # to store patient ID with parotids
for key in structures_name.keys():
    structures_name_clean[key] = 0
    for key2 in structures_name[key].keys():
        if 'parotid' in structures_name[key][key2].lower():
            structures_name_clean[key] = 1
            continue
Struct_data_parotid = pd.DataFrame(structures_name_clean)
Struct_data_parotid.to_csv(os.path.join(path,"structure_names.csv"))

for id in structures_name_clean.keys():
    if structures_name_clean[id] == 0:
        print(id)

        
# flag_parotid = [] # to extract which patient has contoured parotids
# flag_mandible = [] # to extract which patient has contoured mandible
# for key in structures_name.keys():
#     for key2 in structures_name[key].keys():
#         if 'parotid' in structures_name[key][key2].lower():
#             flag_parotid.append(key)
#             print(key)
#         elif 'mandible' in structures_name[key][key2].lower():
#             flag_mandible.append(key)


#####################################################################################
####################### Calculate DVH of all parotids ###############################

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

