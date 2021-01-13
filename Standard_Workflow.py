
## Step 1 ##
# Rename four structure names(parotid l&r ,spinal cord and brain stem) in data folder 
from Anatomic_Features import Standard_Tools
data_folder = "E:/NBIA_HNSCC_DATA/HNSCC_clean_oropharynx_PTVs/"
input_data_path = "E:/demo/NBIA_Head_Neck_DATA/"
ptv_name = "E:/NBIA_HNSCC_DATA/PTVInformationStatistics_Clean.xlsx"
# ptv_name = "C:/Users/xhuae08006/OneDrive - Elekta/Desktop/PTVInformationStatistics_Checking.xlsx"
Y = Standard_Tools(data_folder,input_data_path,ptv_name)
# Y.Rename_STRTNAME_in_rtssdcm()

## Step 2 ## 
# To move all the rtss data and rtdose data to lonely folders
rtss_folder = "E:/NBIA_HNSCC_DATA/HNSCC_rtss_PTVs/"
rtplan_folder = "E:/NBIA_HNSCC_DATA/HNSCC_rtplan_PTVs/"
rtdose_folder = "E:/NBIA_HNSCC_DATA/HNSCC_rtdose_PTVs/"
# Y.Move_RTSSDIOCM_To_Folder(rtss_folder,rtplan_folder,rtdose_folder)
## Step 3 ## 
# Calculate structure volume of ptv total, parotid l&r, spinal cord and brain stem in data folder (Volumes.csv)
# As I add new structure in data folder, so all structure should be re-tried.

######################################################################################################
######################### To calculate the volume of each structure ##################################
######################################################################################################
save = False
Y.Calculate_Clean_Volume(rtss_folder,save)
Y.Volume_clean_data 
## Step 4 ##
# Calculate DTH for all Total ptv with four OARs (e.g. Parotid_L_DTH.csv)
Y.Batch_Deal_with_DTH(rtss_folder,'Parotid_L')
Y.Batch_Deal_with_DTH(rtss_folder,'Parotid_R')
Y.Batch_Deal_with_DTH(rtss_folder,'SpinalCord')
Y.Batch_Deal_with_DTH(rtss_folder,'BrainStem')

## Step 5 ##
# Calculate DVH for all four OARs (e.g. Parotid_L_DVH.csv)
Y.Batch_Deal_with_DVH(rtss_folder,rtdose_folder,'Parotid_L')
Y.Batch_Deal_with_DVH(rtss_folder,rtdose_folder,'Parotid_R')
Y.Batch_Deal_with_DVH(rtss_folder,rtdose_folder,'SpinalCord')
Y.Batch_Deal_with_DVH(rtss_folder,rtdose_folder,'BrainStem')
Y.Batch_Deal_with_DVH(rtss_folder,rtdose_folder,'GTV')

## STL file generation ## (optional)
import os
import pandas as pd
Input_path = "E:\\NBIA_HNSCC_DATA\\HNSCC_rtss_PTVs\\"
CountourSequenceName = "PTV"
Output_path = "E:\\NBIA_HNSCC_DATA\PTV_meshs\\"
ptv_name = pd.read_excel(ptv_name)
kk = ptv_name.iloc[:,1]
kk.index = ptv_name.iloc[:,0]
name_already = [j.split(".")[0] for j in os.listdir(Output_path)]
for item in os.listdir(Input_path):
    if item.split(".")[0] not in name_already and item.split(".")[0] != "HNSCC-01-0047_rtss" and item.split(".")[0] != "HNSCC-01-0057_rtss"\
    and item.split(".")[0] != "HNSCC-01-0067_rtss" and item.split(".")[0] != "HNSCC-01-0072_rtss" and item.split(".")[0] != "HNSCC-01-0074_rtss"\
    and item.split(".")[0] != "HNSCC-01-0090_rtss" and item.split(".")[0] != "HNSCC-01-0091_rtss" and item.split(".")[0] != "HNSCC-01-0097_rtss"\
    and item.split(".")[0] != "HNSCC-01-0134_rtss" and item.split(".")[0] != "HNSCC-01-0176B_rtss":
        print("The item is {}".format(item.split(".")[0]))
        InputFileName = os.path.join(Input_path,item)
        OutputFileName = os.path.join(Output_path,item.split(".")[0]+".stl")
        CountourSequenceName = kk[item.split(".")[0]]
        Y.RTSS_DCM_TO_MESH(InputFileName,CountourSequenceName,OutputFileName)

## display GTV DVH to show the prescription dose.
import pandas as pd

GTV_DVH = pd.read_csv("D:\\demo\\NBIA_Head_Neck_DATA\\DVH_GTV_Clean.csv")

