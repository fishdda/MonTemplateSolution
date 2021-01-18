# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 11:42:10 2020

@author: Henry Huang in Elekta Shanghai Co. Ltd.
"""
       
from MONACO511_TG_WEB import Initialization_MON511
# from GT_MONACO_20200312 import Initialization_MON551
import os
###############################################################################
#################################   GUI   #####################################
pt_id = '002'
#,'003','004','005','006','0010','0016','0020'
delivery_method = 'VMAT'
# pt_id_list = ['002']
# delivery_method_list = ['IMRT','VMAT']
fx = 28
prep_dose=61.6
grid_dose = 2.5
###############################################################################


# default file path 
path = 'C:/Users/xhuae08006/OneDrive - Elekta/Documents/MonTemplateSolution/HYPSolution511'

# absolute path for electronic protocol 
protocol_xlsx = os.path.join(path,'XH protocol.xlsx')

# absolute path for structure name changes
PT_path = 'C:/Users/Public/Documents/CMS/FocalData/Installation/5~Clinic_XH/1~'


X = Initialization_MON511(pt_id,
                delivery_method,
                fx,
                prep_dose,
                grid_dose,
                path,
                protocol_xlsx)
        
X.MAIN_GENERATE('NPC')
###############################################################################