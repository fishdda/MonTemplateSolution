
from DataProcess import DATAPROCESS
import pymedphys
import numpy as np
CT_path = 'C:/GitFolder/RL-Application-in-TPS/AUTO-PLANNING/AutoTemplateTuning/projects/dose prediction/DATA'
DATA_path = 'C:/GitFolder/RL-Application-in-TPS/AUTO-PLANNING/AutoTemplateTuning/projects/sequencing prediction/DATA_/'

X = DATAPROCESS(CT_path,DATA_path)
# PLAN = X.TRANSFER_RAW_DATA()

TRF_path = 'D:\\VM Settings\\XHTOMODATA RAW\\TRF\\005.trf'
DICOM_path = 'D:\\VM Settings\\XHTOMODATA RAW\\005\\005_VMATQATEST.dcm'


DICOM_MU_DENSITY = X.MU_Density_DICOM(DICOM_path) # DICOM path
np.save(DATA_path + '005_MU_density.npy',DICOM_MU_DENSITY)
TRF_MU_DENSIT = X.MU_Density_TRF(TRF_path) # TRF path
np.save(DATA_path + '005_MU_density_trf.npy',TRF_MU_DENSIT)
# grid = pymedphys.mudensity.grid() 
# pymedphys.mudensity.display(grid,TRF_MU_DENSIT) # plot the TRF MU density             ###
# pymedphys.mudensity.display(grid,DICOM_MU_DENSITY) # plot the DICOM MU density             ###

# display all
X.display_mu_density_and_difference()
                                      ###
# calculate gamma distribution
import numpy as np
import os
import matplotlib.pyplot as plt
pt_id = ['002','003','004','005','006','009','0010','0015','0016','0017','0019','0020','0023','0024','0025','0026','0028']
DATA_path = 'C:/GitFolder/RL-Application-in-TPS/AUTO-PLANNING/AutoTemplateTuning/projects/sequencing prediction/DATA_/'
# diff,diff_mean = [],[]
# for item in pt_id:
#     DICOM_MU = np.load(os.path.join(DATA_path,item + '_MU_density.npy'))
#     TRF_MU = np.load(os.path.join(DATA_path,item + '_MU_density_trf.npy'))
#     D = DICOM_MU - TRF_MU
#     diff.append(D.flatten())
#     diff_mean.append(np.mean(np.abs(D)))

# fig = plt.figure(1, figsize=(9, 6))
# ax = fig.add_subplot(111)
# bp = ax.boxplot(diff,patch_artist=True)
# ax.set_xticklabels(['002','003','004','005','006','009','0010','0015','0016','0017','0019','0020','0023','0024','0025','0026','0028'])
# ax.set_xlabel('VMAT case #')
# ax.set_ylabel('MU difference between DICOM and TRF')
# plt.show()

# import pydicom
# plan_path = 'C:\\GitFolder\\RL-Application-in-TPS\\AUTO-PLANNING\\AutoTemplateTuning\\projects\\dose prediction\\DATA\\002\\002_VMATQATEST_Dose.dcm'
# plan = pydicom.read_file(plan_path,force=True)

# strt_path = 'C:\\GitFolder\\RL-Application-in-TPS\\AUTO-PLANNING\\AutoTemplateTuning\\projects\\dose prediction\\DATA\\002\\002_StrctrSets.dcm'
# strt = pydicom.read_file(strt_path,force=True)

# test_dvh = plan.DVHSequence[0].DVHData
# test_dvh = [float(item) for item in test_dvh]

# X = test_dvh[::2]
# Y = test_dvh[1::2]

# DVH curves
