# This Script used for analyze anatomy features to 
# 1) calculate DTH to represent the Volume Overlap between target and OARs
# 2) calculate the volume of each structure.
# 3) Deep
class Anatomy_Features:

    def __init__(self,rtss_path):
        self.rtss_path = rtss_path
        
    def _find_nearest_index(self,array, value):
        ''' Find the nearest index of Z slices.

        Args:
            array: a list of Z positions
            value: a value want to be found
           
        '''
        n = [abs(i-value) for i in array]
        idx = n.index(min(n))
        return idx
    
    def _point_inside_polygon(self,x,y,poly):
        ''' To determine if the point (x,y) in a plane

        Args:
            x,y: the corresponding x position and y position in
                 the plane
            poly: all the point in this planes(list)
          
        '''
        n = len(poly)
        inside =False
        
        if poly == []:
            
            inside =False
            
        else:
            
            p1x,p1y = poly[0]
            for i in range(n+1):
                p2x,p2y = poly[i % n]
                if y > min(p1y,p2y):
                    if y <= max(p1y,p2y):
                        if x <= max(p1x,p2x):
                            if p1y != p2y:
                                xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                            if p1x == p2x or x <= xinters:
                                inside = not inside
                p1x,p1y = p2x,p2y

        return inside

    def _calculate_DTH(self, rtss_path, Target_ID, OAR_ID, wplt=False):
        ''' To calculate the distance to target histogram of OARs

        Args:
            Target_ID : the corresponding PTV or CTV ID stored in rtss.dcm file
            OAR_ID    : the corresponding OARS ID stored in rtss.dcm
            wplt(True/False): whether to plot the DTH of OARs
           
        '''
        from dicompylercore import dicomparser
        import numpy as np
        import matplotlib.pyplot as plt

        rtss = dicomparser.DicomParser(rtss_path)
        structures = rtss.GetStructures()

        ## Target Information (PTV)
        PTV = structures[Target_ID]
        PTV['planes'] = rtss.GetStructureCoordinates(Target_ID)
        key_PTV = PTV['planes'].keys() # slice of PTV
        try:
            z_ptv_min = min([float(j) for j in key_PTV])
            z_ptv_max = max([float(j) for j in key_PTV])
            ptv_points=[]
        except ValueError:
            print("min() arg is an empty sequence")
            return {'DTH':[]}
        for j in key_PTV:
            ptv_points.extend(PTV['planes'][j][0]['data'])
        ptv_points = np.array(ptv_points)

        ## #Non-PTV Information
        OARs = structures[OAR_ID]
        try:
            OARs['planes'] = rtss.GetStructureCoordinates(OAR_ID)
            key_OARs = OARs['planes'].keys() # the slice of OARs

            OARs_points = []
            poly_PTV = {}
            for j in key_OARs:
                Z_OARs = OARs['planes'][j][0]['data'][0][2] # enumerate Z slice
                # Find if the Z slice is within the PTV
                if Z_OARs >= z_ptv_min and Z_OARs <= z_ptv_max:
                    
                    Temp_Data = PTV['planes'][list(key_PTV)[self._find_nearest_index([float(jss) for jss in key_PTV],Z_OARs)]][0]['data']
                    poly_PTV[Z_OARs] = [[Temp_Data[k][0],Temp_Data[k][1]]for k in range(len(Temp_Data))] # the correponding ptv planes
                else:
                    poly_PTV[Z_OARs] = []
            
                ## to search the min x,y positions and max x,y positions in this Z slice
                x_min = OARs['planes'][j][0]['data'][0][0]
                x_max = OARs['planes'][j][0]['data'][0][0]
                y_min = OARs['planes'][j][0]['data'][0][1]
                y_max = OARs['planes'][j][0]['data'][0][1]

                for k in range(0,len(OARs['planes'][j][0]['data'])):
                    if x_min > OARs['planes'][j][0]['data'][k][0]:
                        x_min = OARs['planes'][j][0]['data'][k][0]
                    if x_max < OARs['planes'][j][0]['data'][k][0]:
                        x_max = OARs['planes'][j][0]['data'][k][0]
                    if y_min > OARs['planes'][j][0]['data'][k][1]:
                        y_min = OARs['planes'][j][0]['data'][k][1]
                    if y_max < OARs['planes'][j][0]['data'][k][1]:
                        y_max = OARs['planes'][j][0]['data'][k][1]
                    del OARs['planes'][j][0]['data'][k][2]
                poly_OARs = OARs['planes'][j][0]['data']

                x_min = int(x_min)
                y_min = int(y_min)
                x_max = int(x_max)
                y_max = int(y_max)

                if x_min != x_max and y_min != y_max:
                    for x in range(x_min,x_max):
                        for y in range(y_min,y_max):
                            if self._point_inside_polygon(x,y,poly_OARs):
                                OARs_points.append([x,y,Z_OARs])
            
            OARs_points = np.array(OARs_points)
            S = [np.min(np.sqrt(np.sum((ptv_points - OARs_points[i])*(ptv_points - OARs_points[i]),1))) \
                for i in range(OARs_points.shape[0])]
            
            for kk in range(len(S)):
                if self._point_inside_polygon(OARs_points[kk][0],OARs_points[kk][1],poly_PTV[OARs_points[kk][2]]):
                    S[kk] = -S[kk]
            try:
                X = np.round(np.arange(np.min(S),np.max(S),0.01),2)
                S = np.round(S,2)
                Y = [np.sum(S == i) for i in X]
                Y1 = np.cumsum(Y)/np.sum(Y)               
            except ValueError:
                print("ValueError: zero-size array to reduction operation minimum which has no identity")
                return {'DTH':[]}
            # print('X shapes:{},Y shapes:{}'.format(X.shape,Y.shape))

            if wplt:
                print('Figure!')
                fig = plt.figure()
                axes= fig.add_axes([0.1,0.1,0.8,0.8])
                axes.plot(X,Y1)
                plt.xlabel('distance(mm)')
                plt.ylabel('volume(%)')
                plt.title('DTH_'+structures[OAR_ID]['name']+'_'+structures[Target_ID]['name'])
                axes.set_ylim([0,1])
                plt.grid(b=True)
                plt.show()
                # plt.savefig('DTH_'+structures[OAR_ID]['name']+'.png',fig)

                return {'DTH':[X,Y1]}
            else:
                print('No figure was shown!')
                return {'DTH':[X,Y1]}

        except AttributeError:
            print("'Dataset' object has no attribute 'ContourGeometricType'")
            return {'DTH':[]}

    def _calculate_polygon_area(self,polydata):
        '''
           This function was used to calculate area of polygon
           https://www.wikihow.com/Calculate-the-Area-of-a-Polygon
           
        '''
        area = 0
        for i,item in enumerate(polydata):

            if i+1 >= polydata.shape[0]:
                area += (polydata[i,0]*polydata[0,1]) - (polydata[0,0]*polydata[i,1])
            else:
                area += (polydata[i,0]*polydata[i+1,1]) - (polydata[i+1,0]*polydata[i,1])

        area = area*0.5

        return abs(area)

    def _calculate_structure_volume(self,ptv_points):
        '''
           This function was used to calculate the volume of a structure 
           (To calculate the contour area is a smart algorithm)

        Args: 
            ptv_points: the contour points in PTV or any structure stored in rtss.dcm
            
        '''
        import numpy as np
        volume = 0
        slice = np.unique(ptv_points[:,-1])
        slice_thickness = slice[1] - slice[0]
        for i,item in enumerate(slice):

            slice_data = ptv_points[ptv_points[:,-1] == item]

            volume += self._calculate_polygon_area(slice_data[:,0:-1])*slice_thickness
            #print('slice:{},area:{}'.format(item,self._calculate_polygon_area(slice_data[:,0:-1])))

        return volume 

    def DVH_data_extration(self,data_path,pt_id,structure_name):
        '''
          This function was used for DVH data extraction for display and
          further analysis
        '''

        # test_path = 'C:\\GitFolder\\RL-Application-in-TPS\\AUTO-PLANNING\\AutoTemplateTuning\\projects\\DVH_prediction\\DATA_'
        import numpy as np
        import pandas as pd
        import os
        ParotidR = {}
        length = []
        for item in pt_id:
            Temp_DVH = np.load(os.path.join(data_path,item+'_DVHDATA.npy'),allow_pickle=True)
            Temp_DVH = Temp_DVH.item()
            Parotid_R = Temp_DVH[structure_name]
            Parotid_R = [float(item) for item in Parotid_R]
            X = [item for i,item in enumerate(Parotid_R) if i%2 == 0]
            X = np.cumsum(np.array(X))
            ParotidR[item] = [item for i,item in enumerate(Parotid_R) if i%2 == 1]
            length.append(len(ParotidR[item]))
            
        SS = max(length)
        DATA = np.zeros([len(pt_id),SS])
        for i,item in enumerate(pt_id):
            DATA[i,0:len(ParotidR[item])] = ParotidR[item]

        DATA = pd.DataFrame(DATA.T,columns = pt_id)
        new_X = list(2*np.arange(DATA.index[0],DATA.index[-1]+1))
        DATA.index = new_X
        DATA_new = DATA.T

        return DATA_new

    def DTH_data_extration(self,DTH_STAT_DATA,pt_id):
        '''
           This function was prepared for DTH data extraction and
           preparation for later usage
        '''
        from scipy import interp
        import numpy as np
        import pandas as pd
        X = [(np.min(DTH_STAT_DATA[i][0]),np.max(DTH_STAT_DATA[i][0])) for i in range(len(DTH_STAT_DATA))]
        X_min = np.min(np.array([item[0] for item in X]))
        X_max = np.max(np.array([item[1] for item in X]))
        X_sum = np.arange(X_min,X_max,0.01)
        y1_new = []
        for i in range(len(DTH_STAT_DATA)):
            y1_new.append(interp(X_sum,DTH_STAT_DATA[i][0],DTH_STAT_DATA[i][1]))
        y1_new = np.array(y1_new)
        y_new = pd.DataFrame(y1_new.T,index = X_sum,columns = pt_id)

        return y_new





class Standard_Tools:
    def __init__(self,data_folder,input_data_path,ptv_std_name):
        self.data_folder = data_folder
        self.written_path = input_data_path
        self.ptv_std_name = ptv_std_name

    def Rename_STRTNAME_in_rtssdcm(self):
        '''
          Rename Parotid Left, Right, SpinalCord and BrainStem in data
          DICOM DATA:
          -- CT1.dcm
          -- CT2.dcm
          ...
          -- CTn.dcm
          -- rtss.dcm  (Modified Names in rtss.dcm)
          -- rtdose.dcm
          -- rtplan.dcm

        '''
        import pydicom
        import os
        # input_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_rtss\\"
        # output_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_rtss\\"
        # dirty_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_rtss\\"
        clean_rtss_name = os.listdir(self.data_folder)
        for name in clean_rtss_name:
            print("Name:{}".format(name))
            dcm = pydicom.read_file(os.path.join(self.data_folder,name,name+'_rtss.dcm'),force=True)
            count = 0
            for ROI in dcm.StructureSetROISequence:
                ## Brain Stem and Spinal Cord parts ##
                if 'expanded' not in ROI.ROIName.lower() and \
                    'expd' not in ROI.ROIName.lower() and \
                        'exp' not in ROI.ROIName.lower() and \
                        'hot' not in ROI.ROIName.lower() and \
                    'expand' not in ROI.ROIName.lower() and \
                        'new' not in ROI.ROIName.lower() and \
                        'prv' not in ROI.ROIName.lower() and \
                    '0.5cm' not in ROI.ROIName.lower() and \
                    'cords' not in ROI.ROIName.lower() and \
                        'no' not in ROI.ROIName.lower() and \
                        'mp' not in ROI.ROIName.lower() and \
                        'pcv' not in ROI.ROIName.lower() and \
                    'avoid' not in ROI.ROIName.lower() and \
                    '0.5mm' not in ROI.ROIName.lower() and \
                'planning' not in ROI.ROIName.lower() and \
                ' 2' not in ROI.ROIName.lower():
                    if 'cord' in ROI.ROIName.lower():
                        print('original name:{}'.format(ROI.ROIName)) 
                        dcm.StructureSetROISequence[dcm.StructureSetROISequence.index(ROI)].ROIName = 'SpinalCord'
                        print('new name:{}'.format(dcm.StructureSetROISequence[dcm.StructureSetROISequence.index(ROI)].ROIName))
                        count +=1
                    if 'stem' in ROI.ROIName.lower():
                        print('original name:{}'.format(ROI.ROIName))
                        dcm.StructureSetROISequence[dcm.StructureSetROISequence.index(ROI)].ROIName = 'BrainStem'
                        print('new name:{}'.format(dcm.StructureSetROISequence[dcm.StructureSetROISequence.index(ROI)].ROIName))
                        count +=1

                ## Parotid left and Parotid right parts ##
                if 'parotid' in ROI.ROIName.lower() and 'sub' not in ROI.ROIName.lower() and \
                'avoid' not in ROI.ROIName.lower() and 'low' not in ROI.ROIName.lower() and \
                '20' not in ROI.ROIName.lower() and 'tail' not in ROI.ROIName.lower() and \
                'sud' not in ROI.ROIName.lower() and 'ptv' not in ROI.ROIName.lower() and \
                'parotid2' not in ROI.ROIName.lower() and 'parotids' not in ROI.ROIName.lower() and \
                    'push' not in ROI.ROIName.lower() and 'hot' not in ROI.ROIName:
                    if 'lt' in ROI.ROIName.lower() or 'left' in ROI.ROIName.lower() or 'L' in ROI.ROIName or 'l ' in ROI.ROIName or 'lft' in ROI.ROIName:
                        print('original name:{}'.format(ROI.ROIName))
                        dcm.StructureSetROISequence[dcm.StructureSetROISequence.index(ROI)].ROIName = 'Parotid_L'
                        print('new name:{}'.format(dcm.StructureSetROISequence[dcm.StructureSetROISequence.index(ROI)].ROIName))
                    elif 'rt' in ROI.ROIName.lower() or 'right' in ROI.ROIName.lower() or 'R' in ROI.ROIName or 'r ' in ROI.ROIName:
                        print('original name:{}'.format(ROI.ROIName))
                        dcm.StructureSetROISequence[dcm.StructureSetROISequence.index(ROI)].ROIName = 'Parotid_R'
                        print('new name:{}'.format(dcm.StructureSetROISequence[dcm.StructureSetROISequence.index(ROI)].ROIName))

            print("Number:{}".format(count))
            dcm.save_as(os.path.join(self.data_folder,name,name+"_rtss.dcm"))

    def Move_RTSSDIOCM_To_Folder(self,rtss_folder,rtplan_folder,rtdose_folder):
        '''
          To Move all rtss dicom to a specific folder to store all the rtss.dcm
          to make further analysis easier and more efficient
        '''
        import os
        import pydicom
        #################################################################################################
        ############################ To put all the rtss.dcm into a folder ##############################
        #################################################################################################
        # data_folder = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_PTVs\\"
        # HNSCC_clean_oropharynx_PTVs
        # -- HNSCC-01-0001
        # ----- CT1.dcm
        # ----- CT2.dcm
        # ----- rtss.dcm
        # ----- rtdose.dcm
        # ----- rtplan.dcm
        data_name = os.listdir(self.data_folder)
        for i,item in enumerate(data_name):
            print("HNSCC number:{} and name: {}".format(i+1,item))
            sub_folder = os.listdir(os.path.join(self.data_folder,item))
            # sub_folder_names = os.listdir(os.path.join(self.data_folder,item,sub_folder[0]))
            # for j,jtem in enumerate(sub_folder):
            #     sub_sub_folder_names = os.listdir(os.path.join(self.data_folder,item,sub_folder,jtem))
            for file_name in sub_folder:
                dcm = pydicom.read_file(os.path.join(self.data_folder,item,file_name),force=True)
                if dcm.Modality == 'RTSTRUCT':
                    print("DCM Modality:{}".format(dcm.Modality))
                    dcm.save_as(os.path.join(rtss_folder,item,item+'_rtss.dcm'))
                elif dcm.Modality == 'RTDOSE':
                    print("DCM Modality:{}".format(dcm.Modality))
                    dcm.save_as(os.path.join(rtdose_folder,item,item+'_rtdose.dcm'))
                elif dcm.Modality == 'RTPLAN':
                    print("DCM Modality:{}".format(dcm.Modality))
                    dcm.save_as(os.path.join(rtplan_folder,item,item+'_rtplan.dcm'))
                elif dcm.Modality == 'CT':
                    continue       

    def Calculate_Clean_Volume(self,rtss_folder,save):
        '''
        Calculate the Clean Volume of each structures like PTV, Parotid L&R, Spinal Cord, Brain Stem.
        '''
        import os 
        import math
        import pydicom
        import numpy as np
        import pandas as pd
        from Anatomic_Features import Anatomy_Features
        X = Anatomy_Features(rtss_folder)
        file_names = os.listdir(rtss_folder)
        Volume = {item.split('.')[0]:{} for item in file_names}
        # key = "HNSCC-01-0004_rtss"
        ptv_name = pd.read_excel(self.ptv_std_name)
        kk = ptv_name.iloc[:,1]
        kk.index = ptv_name.iloc[:,0]  # HNSCC-01-0004_rtss:Total

        self.NaN_List = []

        for key in kk.index:
            rtss = pydicom.read_file(os.path.join(rtss_folder,key+".dcm"),force=True)
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
                if structures_name[Contour.ReferencedROINumber] == kk[key]:
                    if math.isnan(vol_cc):
                        self.NaN_List.append((key,vol_cc))
                    else:
                        Volume[key]['PTV_Combinations'] = vol_cc
                else:
                    Volume[key][structures_name[Contour.ReferencedROINumber]] = vol_cc
            
            print("pt:{}, and num of structures:{}".format(key,len(Volume[key].keys())))

        self.Volume_data = pd.DataFrame(Volume)
        self.Volume_data = self.Volume_data.T
        self.Volume_clean_data = self.Volume_data[['Parotid_L','Parotid_R','SpinalCord','BrainStem','PTV_Combinations']]
        if save == True:
            self.Volume_clean_data.to_csv(os.path.join(rtss_folder,"Volume_Clean_New.csv"))


    def Batch_Deal_with_DTH(self,rtss_path,strt_name):
        '''
           This script is used for dealing with DTH data including 
           cleaning, renaming and write to csv
           strt_name: the name of structure
        '''
        DTH_Left,DTH_Right = {},{}
        left_log_inf = []
        import os
        import numpy as np
        import pandas as pd
        from Anatomic_Features import Anatomy_Features
        from dicompylercore import dicomparser,dvh, dvhcalc
        X = Anatomy_Features(1)
        # rtss_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_rtss\\"
        # dose_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_rtdose\\"
        # log_path = "C:\\Users\\xhuae08006\\OneDrive - Elekta\\Desktop\\demo\\"
        ptv_name = pd.read_excel(self.ptv_std_name)
        kk = ptv_name.iloc[:,1]
        kk.index = ptv_name.iloc[:,0]
        filenames = os.listdir(rtss_path)
        for name in filenames:
            ## calculate Parotid DTH
            DTH_Left[name.split('.')[0]] = {}
            DTH_Right[name.split('.')[0]] = {}
            rtss = dicomparser.DicomParser(os.path.join(rtss_path,name))
            rtstructures = rtss.GetStructures()
            Target_ID,OAR_ID = {},{}
            for i in rtstructures.keys():
                if rtstructures[i]['name'].lower() == kk[name.split('.')[0]].lower():
                    Target_ID[rtstructures[i]['name']] = rtstructures[i]['id']
                    print("xlsx name:{} and rtstructures[i]['name']:{}".format(kk[name.split('.')[0]],rtstructures[i]['name']))
                elif strt_name in rtstructures[i]['name']:
                    OAR_ID[rtstructures[i]['name']] = rtstructures[i]['id']
                    
            for key1 in Target_ID.keys():
                for key2 in OAR_ID.keys():
                    if key2 == strt_name:
                        print("PT {} for DTH is calculated between {} and {}".format(name.split('.')[0],key1,key2))
                        DTH_data = X._calculate_DTH(os.path.join(rtss_path,name),Target_ID[key1], OAR_ID[key2], wplt=False)
                        DTH_Left[name.split('.')[0]][key1] = DTH_data
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
        DTH_Left_Matrix.to_csv(os.path.join(self.written_path,"DTH_"+strt_name+"_Clean.csv"))


    def Batch_Deal_with_DVH(self,rtss_path,dose_path,strt_name):
        '''
           This script is used for dealing with DVH data including cleaning, renaming and write to csv
           strt_name: the name of structure
        '''
        ## calculate DVH ##
        ## Just try to calculate DTH all the cases
        DVH = {}
        log_inf = []
        import os
        import pandas as pd
        import numpy as np
        from Anatomic_Features import Anatomy_Features
        from dicompylercore import dicomparser,dvh, dvhcalc
        # X = Anatomy_Features(1)
        # rtss_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_rtss\\"
        # dose_path = "D:\\NBIA_HNSCC_DATA\\HNSCC_clean_oropharynx_rtdose\\"
        # log_path = "C:\\Users\\xhuae08006\\OneDrive - Elekta\\Desktop\\demo\\"
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
                elif strt_name in rtstructures[i]['name']:
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
                if  key2 == strt_name:
                    DVH_Left_Parotid[key] = DVH[key][key2]



        Left_Parotid_Dose_max =[(DVH_Left_Parotid[key]['Dose'][-1],key) for key in DVH_Left_Parotid.keys()]
        # Right_Parotid_Dose_max = [(DVH_Right_Parotid[key]['Dose'][-1],key) for key in DVH_Right_Parotid.keys()]
        Left_Parotid_Dose_Index = np.arange(0,max(Left_Parotid_Dose_max)[0],0.01)
        # Right_Parotid_Dose_Index = np.arange(0,max(Right_Parotid_Dose_max)[0],0.01)

        Left_Parotid_DVH_Matrix = np.zeros([len(DVH_Left_Parotid.keys()),Left_Parotid_Dose_Index.shape[0]]) 
        # Right_Parotid_DVH_Matrix = np.zeros([len(DVH_Right_Parotid.keys()),Right_Parotid_Dose_Index.shape[0]]) 
        columns_left,columns_right = [],[]
        for i,item in enumerate(DVH_Left_Parotid.keys()):
            Left_Parotid_DVH_Matrix[i,0:DVH_Left_Parotid[item]['Volume'].shape[0]] = DVH_Left_Parotid[item]['Volume']
            columns_left.append(item)
        # for i,item in enumerate(DVH_Right_Parotid.keys()):
            # Right_Parotid_DVH_Matrix[i,0:DVH_Right_Parotid[item]['Volume'].shape[0]] = DVH_Right_Parotid[item]['Volume']
            # columns_right.append(item)

        
        Left_Parotid_DVH_Matrix_ = pd.DataFrame(Left_Parotid_DVH_Matrix.T,columns=columns_left,index = Left_Parotid_Dose_Index)
        # Right_Parotid_DVH_Matrix_ = pd.DataFrame(Right_Parotid_DVH_Matrix.T,columns=columns_right,index = Right_Parotid_Dose_Index)
        Left_Parotid_DVH_Matrix_.to_csv(os.path.join(self.written_path,"DVH_"+strt_name+"_Clean.csv"))
        # Right_Parotid_DVH_Matrix_.to_csv("C:\\Users\\xhuae08006\\OneDrive - Elekta\\Desktop\\demo\\NBIA_Head_Neck_DATA\\DVH_SpinalCord_Clean.csv")

    def RTSS_DCM_TO_MESH(self,InputFileName,CountourSequenceName,OutputFileName):
        '''
          This Script is used for generating mesh data from rtss dicom files
          it was copied from 2016-2018 Lukasz J. Nowak

          InputFileName: the rtss.dcm file path used for extracting contour and structures (001_rtss.dcm)
          CountourSequenceName: specific structure name ("PTV")
          OutputFileName: the output folder path for saving the rtss.dcm
        '''
        import os
        #start with importing pydicom library - it will be used for reading specific content from the DICOM file
        import pydicom as dicom
        # path = 'D:/NBIA_HNSCC_DATA/HNSCC_clean_oropharynx_PTVs/'
        # pt_id = 'HNSCC-01-0168'
        # InputFileName = os.path.join(path,pt_id,pt_id+'_rtss.dcm')
        # CountourSequenceName='Total ptv'
        # OutputFileName = 'D:/'+pt_id+'_'+CountourSequenceName+'.stl'

        
        # InputFileName = "C:\\Users\\Public\\Documents\\CMS\\FocalData\\DCMXprtFile\\111111_StrctrSets.dcm"
        # CountourSequenceName = "PTV"
        # OutputFileName = 'D:/prostate.stl'




        #Class definitions:
        #*********************************
        #Class of points:
        class point:
            def __init__(self,x,y,z):
                self.x=x
                self.y=y
                self.z=z
            #Two points are considered equal, if all their coordinates are equal:
            def __eq__(self,other):
                if self.x == other.x and self.y == other.y and self.z == other.z:
                    return True
                else:
                    return False


        #*********************************
        #stlfacet class defines structures acordingly to the STL format. 
        #Triangle elements with vertices numbered in such a way, that the normal versor points outwards    
        class stlfacet:
            def __init__(self,pointa,pointb,pointc):
                self.pointa = pointa
                self.pointb = pointb
                self.pointc = pointc     
                #coefficients of vector normal to the facet surface are computed using point class:
                vectu = point(pointb.x-pointa.x,pointb.y-pointa.y,pointb.z-pointa.z)
                vectv = point(pointc.x-pointa.x,pointc.y-pointa.y,pointc.z-pointa.z)
                nnx = vectu.y*vectv.z - vectu.z*vectv.y
                nny = vectu.z*vectv.x - vectu.x*vectv.z
                nnz = vectu.x*vectv.y - vectu.y*vectv.x
                #normalization:
                self.nx = nnx / ((nnx**2 + nny**2 + nnz**2)**(1/2))
                self.ny = nny / ((nnx**2 + nny**2 + nnz**2)**(1/2))
                self.nz = nnz / ((nnx**2 + nny**2 + nnz**2)**(1/2))

            #return text with computed coordinates accordingly to the STL file standard:    
            def printfacet(self):
                return 'facet normal ' + str(self.nx) + ' ' + str(self.ny) + ' ' + str(self.nz) + '\n\t outer loop\n\t\t vertex ' + str(self.pointa.x) + ' ' + str(self.pointa.y) + ' ' + str(self.pointa.z) + '\n\t\t vertex ' + str(self.pointb.x) + ' ' + str(self.pointb.y) + ' ' + str(self.pointb.z) + '\n\t\t vertex ' + str(self.pointc.x) + ' ' + str(self.pointc.y) + ' ' + str(self.pointc.z) + '\n\t endloop\n endfacet\n'



        #Function definitions:
        #2D distance between two points:
        def distance2d(point1,point2):
            return ((point2.x-point1.x)**2 + (point2.y-point1.y)**2) ** (1/2)


        #*******************
        #3D distance between two points:
        def distance(point1,point2):
            return ((point2.x-point1.x)**2 + (point2.y-point1.y)**2 + (point2.z-point1.z)**2) ** (1/2)


        #**************************
        #Function findDirection, defined for three subsequent points at a given curve (with middle point having extreme coordinates among all the points at the curve) returns either True or False, depending on the curve orientation
        #True - direction "left", false - direction "right" (accordingly to numeration of points)
        def findDirection(pointa, pointb, pointc): 
            if ((pointb.x - pointa.x) * (pointc.y-pointa.y) - (pointc.x - pointa.x) * (pointb.y - pointa.y)) != 0:
                return ((pointb.x - pointa.x) * (pointc.y-pointa.y) - (pointc.x - pointa.x) * (pointb.y - pointa.y)) < 0
            else:
                return None   #if points are collinear, direction cannot be determined -> value None is returned


        #*******************************
        #Write data to the STL file (input: filename and list of facets)
        def printstl(InputFileName,facets):
            file = open(InputFileName,'w')
            file.write('solid DICOM_contour_model\n')
            for licz in range(len(facets)):
                file.write(facets[licz].printfacet())
            file.write('endsolid DICOM_contour_model\n')
            file.close()


        ##STEP 1
        #Read DICOM file:
        ds=dicom.read_file(InputFileName,force=True)

        #STEP 2: find the specified structure in the file (ROI contour sequence with a given name):
        if not hasattr(ds, 'StructureSetROISequence'):
            print('The specified file does not contain contour sequence data.')
            raise NameError('There is no contour sequence data within the specified DICOM file!')
        for licz in range(len(ds.StructureSetROISequence)):
            if ds.StructureSetROISequence[licz][0x3006,0x26].value==CountourSequenceName:
                contno=licz
                break
            elif licz==len(ds.StructureSetROISequence)-1:
                print('No such contour. Available contour names:')
                for liczc in range(len(ds.StructureSetROISequence)):
                    print(ds.StructureSetROISequence[liczc][0x3006,0x26])
                raise NameError('There is no contour with the specified name within the specified DICOM file!')

        #STEP 3: determine how many slices (curves) are within the specified sequence:
        howManySlices=len(list(ds[0x3006,0x39][contno][0x3006,0x40]))

        #STEP 4: determine number of points in each slice (curve):
        howManyPoints=[0]*howManySlices
        for licz in range(howManySlices):
            howManyPoints[licz]=len(list(ds[0x3006,0x39][contno][0x3006,0x40][licz][0x3006,0x50]))//3  #divide by 3, because each point has 3 coordinates    

        #STEP 5:
        #Read and save the coordinates of each point, in each slice (curve):
        points=[0]*howManySlices   #initialization, first dimension (slices)
        for licz in range(howManySlices):
            points[licz]=[0]*howManyPoints[licz]   #initialization, second dimension [slices][points]
        
        #save values to the initialized lists:  
        for licz1 in range(howManySlices):
            for licz2 in range(howManyPoints[licz1]):
                points[licz1][licz2]=point(float(ds[0x3006,0x39][contno][0x3006,0x40][licz1][0x3006,0x50][3*licz2]),float(ds[0x3006,0x39][contno][0x3006,0x40][licz1][0x3006,0x50][3*licz2+1]),float(ds[0x3006,0x39][contno][0x3006,0x40][licz1][0x3006,0x50][3*licz2+2]))
        #Sort slices accordingle to the increasing z coordinate (slices are within XY plane):
        points.sort(key=lambda x: x[0].z)
        #Update point count:
        howManyPoints=[0]*howManySlices
        for licz in range(len(points)):
            howManyPoints[licz]=len(points[licz])  


        #STEP 6:
        #Remove redundant slices.
        #This version of software assumes, that only one curve per slice (i.e. with one, specified z value for all points belonging to the curve) is permitted. 
        #If greater number of curves with the same z coordinate are detected, only one curve, with greatest number of points is selected for further processing.
        #Other curves are considered as redundant, and removed. In fact, many of such unwanted artifact appears as a result of automated segmentation and curve drawing processes.
        licz=0
        while licz < len(points)-1:
            if points[licz][0].z == points[licz+1][0].z:
                if howManyPoints[licz] >= howManyPoints[licz+1]:
                    del(points[licz+1])
                    del(howManyPoints[licz+1])
                else:
                    del(points[licz])
                    del(howManyPoints[licz])
                licz = -1
            licz += 1
        howManySlices = len(points) #Number of remaining slices is updated


        #STEP 7:
        #numbers of points in each slice (curve) are re-arranged in such a way, that the starting points (indices 0) are the closes points between the subsequent slices (curves)
        #re-numeration starts at point 0, slice 0 (with the lowest z coordinate)
        pt0indices=[0]*howManySlices

        for licz1 in range(howManySlices-1):
            minPointDistance=distance(points[licz1][0],points[licz1+1][0])   #temporal variable for storing distances
            for licz2 in range(howManyPoints[licz1+1]):
                if distance(points[licz1+1][licz2],points[licz1][0]) < minPointDistance:
                    minPointDistance=distance(points[licz1+1][licz2],points[licz1][0])  #current lowest distance value - update
                    pt0indices[licz1+1]=licz2  
                #if starting point in the slice has different index than 0, numbering is re-arranged:
            if not pt0indices[licz1+1]==0:
                points[licz1+1]=points[licz1+1][pt0indices[licz1+1]:] + points[licz1+1][:pt0indices[licz1+1]]

        #STEP 8
        #direction (i.e. points numeration relative to interior/exterior of the closed curve) of each curve is checked.  If differences between slices are found, numbering of points in non-matching curves is reversed 
        #First, an extreme point is found (the point with maximum x or xy coordinates in the curve):
        max_xy=[0]*howManySlices
        max_yi=[0]*howManySlices
        maxIndex=[0]*howManySlices    
        #point with maximum x coordinate is found:
        for licz1 in range(howManySlices): 
            max_xy[licz1]=points[licz1][0].x    
            max_yi[licz1]=points[licz1][0].y    
            for licz2 in range(1,howManyPoints[licz1]):
                if points[licz1][licz2].x > max_xy[licz1]: 
                    max_xy[licz1]=points[licz1][licz2].x
                    max_yi[licz1]=points[licz1][licz2].y
                    maxIndex[licz1]=licz2  #index of point with maximum x coordinate is stored
                elif points[licz1][licz2].x==max_xy[licz1] and points[licz1][licz2].y > max_yi[licz1]:  #if more points have the same max x coordinate, the point with the highest y coordinate is selected among them:
                    max_xy[licz1]=points[licz1][licz2].x
                    max_yi[licz1]=points[licz1][licz2].y
                    maxIndex[licz1]=licz2  #index of the selected point is stored

        #Now, having the extreme point, direction of each curve can be determined ("left" or "right" - directions of all curves should match each other). Three subsequent points are required (i.e. including neighbours of the extreme point):
        sliceDirections=[0]*howManySlices
        for licz in range(howManySlices):
            if maxIndex[licz]==0:
                pointa=points[licz][howManyPoints[licz]-1] #if extreme point hax index 0, then one of its neighbours is the point with maximum index value inside the curve
            else:
                pointa=points[licz][maxIndex[licz]-1]
            pointb=points[licz][maxIndex[licz]]
            if maxIndex[licz]==howManyPoints[licz]-1:
                pointc=points[licz][0] #if extreme point hax index of maximum value, then one of its neighbours is the point with index 0
            else:
                pointc=points[licz][maxIndex[licz]+1]
            sliceDirections[licz]=findDirection(pointa,pointb,pointc)

        #Next, the direction of each curve with relation to the first (bottommost) slice are checked. If differences are found - the numeration of points in non-matching curves is reversed
        direction=sliceDirections[0]
        for licz in range(howManySlices):
            if sliceDirections[licz] != direction:
                points[licz].reverse()  #after reversing, indices must be shifted, such that 0 will be 0 again, not max index value
                points[licz]=points[licz][howManyPoints[licz]-1:]+points[licz][:howManyPoints[licz]-1]       
                
        #STEP 9
        #Discretization of the lateral surface of the considered 3D structure:
        sidefacets=[]  #list of the elements of the discretized lateral surface

        for licz in range(howManySlices-1):    
            pl1ind=0   #index of point within the lower slice
            pl2ind=0    #index of point within the upper slice
            pl1indnext=1
            pl2indnext=1
            end1 = pl1indnext >= howManyPoints[licz]    #just a precausion in case if any of the slices would consist of one point only
            end2 = pl2indnext >= howManyPoints[licz+1]

            while (not end1) or (not end2):     #...until all points on both curves will be included:

                condition = (distance(points[licz][pl1ind],points[licz+1][pl2indnext]) <= distance(points[licz+1][pl2ind],points[licz][pl1indnext]))    #if the current point from the lower slice is closer to the next point at the upper layer, than vice-versa:
                
                if not end2:   #only if there are still more points to link within the upper layer:  
                    if condition or end1:
                        if not points[licz+1][pl2indnext] == points[licz+1][pl2ind]:    #just a precausion in case if points would be doubled (i.e. multiple points with identical coordinates)
                            sidefacets.append(stlfacet(points[licz][pl1ind],points[licz+1][pl2indnext],points[licz+1][pl2ind])) #curve orientation!
                            pl2ind = pl2indnext
                            if pl2ind == 0: #if the we have came back to the first point of the upper layer:
                                end2 = True
                        pl2indnext+=1
                        if  pl2indnext == howManyPoints[licz+1]:   #if the current index is greater than the maximum index, nex iteration will finish in 0
                            pl2indnext = 0

                if not end1: #until the whole lower layer will be covered:    
                    if (not condition) or end2:
                        if not points[licz][pl1indnext] == points[licz][pl1ind]:    #just a precausion in case if points would be doubled (i.e. multiple points with identical coordinates)    
                            sidefacets.append(stlfacet(points[licz+1][pl2ind],points[licz][pl1ind],points[licz][pl1indnext])) #curve orientation!
                            pl1ind = pl1indnext
                            if pl1ind == 0:
                                end1 = True
                        pl1indnext+=1
                        if pl1indnext == howManyPoints[licz]:
                            pl1indnext = 0

        #STEP 10:
        #Discretization of lower and upper base surfaces of the considered 3D structure:
        direction = sliceDirections[0]  #curve orientation as set previously                  
        lowerCap = []  #memory allocation for lower base
        kd = points[0]  

        #duplicated points are removed:
        for pointkd in kd:
            howManykd = kd.count(pointkd)
            if howManykd > 1:
                for licz in range(howManykd - 1):
                    kd.remove(pointkd)

        licz = 0    #counter of subsequent points within the curve
        kdmark = True   #flag preventing curve twisting
        directiond = direction
        while len(kd) >= 3:     #until at least 3 points are left within the curve...
            if findDirection(kd[licz],kd[licz+1],kd[licz+2]) == directiond:   #if we are at the convex part of curve, we create facet and remove the middle point from further considerations
                if not directiond:    #normal versor to the lower base surface must point downwards
                    lowerCap.append(stlfacet(kd[licz],kd[licz+2],kd[licz+1]))    
                else:
                    lowerCap.append(stlfacet(kd[licz],kd[licz+1],kd[licz+2]))   
                kd.remove(kd[licz+1])
                kdmark = True
            licz += 1  
            if licz >= len(kd) - 2:
                if not kdmark:
                    directiond = not directiond 
                licz = 0    
                kdmark = False  

        #Next, we perform analogous processing for the upper base surface:
        upperCap = []
        kg = points[howManySlices-1]

        ###duplicated points are removed:
        for pointkg in kg:
            howManykg = kg.count(pointkg)
            if howManykg > 1:
                for licz in range(howManykg - 1):
                    kg.remove(pointkg)
                    print('point usuniety z kapsla gornego!')

        licz = 0    
        kgmark = True
        directiong = direction
        while len(kg) >= 3:     #until at least 3 points are left within the curve...
            if findDirection(kg[licz],kg[licz+1],kg[licz+2]) == directiong:   #if we are at the convex part of curve, we create facet and remove the middle point from further considerations
                if not directiong:
                    upperCap.append(stlfacet(kg[licz],kg[licz+1],kg[licz+2]))    
                else:
                    upperCap.append(stlfacet(kg[licz],kg[licz+2],kg[licz+1]))    
                kg.remove(kg[licz+1])
                kgmark  = True
            licz += 1   
            if licz >= len(kg) - 2:
                if not kgmark:
                    directiong = not directiong  
                licz = 0    
                kgmark = False  

        #STEP 11:  
        #We save the discretized structure into the STL file with the specified name, using the defined printstl function:
        printstl(OutputFileName,lowerCap + sidefacets + upperCap)

class PCA_SVR_Model:
    '''
       This class was followed by Zhu et al to predict the DVH curves via a PCA+SVR approach
       which is a machine learning method
    '''
    def __init__(self,DVH_DTH_data_path,strt_name):
        '''
           Args:
           DVH_DTH_data_path: this is a pre-defined folder path for storing DVH and DTH data
           strt_name: this is the structure name user wanna predict its DVH

        '''
        
        self.DVH_DTH_data_path = DVH_DTH_data_path
        self.strt_name = strt_name

    def further_clean_up_for_DTH(self,wplt=False):
        '''
        Argvs:
            1. To further clean up the DTH data for model build. e.g. Parotid L may overlap with more than
            one PTVs and choose which should be selected here.
            2. Output the DVH and DTH plot to corresponding folders.
            
        '''
        import matplotlib.pyplot as plt
        import os
        import numpy as np
        DTH_new_columns = []
        
        for name_DVH in self.DVH_Data.columns:
 
            name_DTH = [item for item in self.DTH_Data.columns if name_DVH in item]
            if wplt == True:
                new_path = os.path.join(self.DVH_DTH_data_path,"Fig",name_DVH+"_Fig")
#                os.mkdir(new_path)    
                #DVH show
                print('Figure!')
                fig = plt.figure()
                # axes= fig.add_axes([0.1,0.1,0.8,0.8])
                self.DVH_Data[name_DVH].plot()
                plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
                plt.xlabel('Dose(Gy)')
                plt.ylabel('volume(%)')
                plt.title(self.strt_name+"_DVH")
                # axes.set_ylim([0,1])
                # plt.xlim(-11.61,35.34)
                plt.grid(b=True)
                fig_name1 = os.path.join(new_path,name_DVH+'_'+self.strt_name+"_DVH.tif")
                plt.savefig(fig_name1,dpi=200, bbox_inches='tight')
                plt.show()
                
                # DTH show
                print('Figure!')
                fig = plt.figure()
                # axes= fig.add_axes([0.1,0.1,0.8,0.8])
                self.DTH_Data[name_DTH].plot()
                plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
                plt.xlabel('distance(mm)')
                plt.ylabel('volume(%)')
                plt.title(self.strt_name+"_DTH")
                # axes.set_ylim([0,1])
                # plt.xlim(-11.61,35.34)
                plt.grid(b=True)
                fig_name = os.path.join(new_path,name_DVH+'_'+self.strt_name+"_DTH.tif")
                plt.savefig(fig_name,dpi=200, bbox_inches='tight')
                plt.show()
             
            ## further extract DTH to make it the same dimension with DVH matrix in patient number
            SS = np.abs(self.DTH_Data[name_DTH].index - 0)
            indx = np.where(SS == np.min(SS))
            Temp = self.DTH_Data[name_DTH].iloc[indx[0][0],:]
            indx_1 = np.where(Temp == np.max(Temp))
            DTH_new_columns.append(Temp.index[indx_1[0][0]])
        
        return DTH_new_columns


    def PCA_Compression(self,num_components,wplt,DVH_file_name,DTH_file_name):

        '''
           This script was used for principle component analysis to 
           compress the data dimension of DTH and DVH to build up the correlationship
           
           Args:
            1) num_components: the principle components selected
            2) wplt: True means plot the principle components and output to a fixed folder path
                     False means no plot and no output
            3) DVH_file_name: DVH data file name selected to compress
            4) DTH_file_name: DTH data file name selected to compress
            5) pt_id: all the NPC cases IDs
        '''
        import os
        import pandas as pd
        from sklearn.preprocessing import StandardScaler
        from sklearn.decomposition import PCA

        ## read the DVH 
        self.DVH_path = os.path.join(self.DVH_DTH_data_path,DVH_file_name)
        self.DTH_path = os.path.join(self.DVH_DTH_data_path,DTH_file_name)

        DVH_Data = pd.read_csv(self.DVH_path)
        DTH_Data = pd.read_csv(self.DTH_path)
        
        self.DVH_Data = DVH_Data.iloc[:,1:]
        self.DTH_Data = DTH_Data.iloc[:,1:]
        
        self.DVH_Data.index = DVH_Data.iloc[:,0]
        self.DTH_Data.index = DTH_Data.iloc[:,0]
        
        ## To drop out the NaN value from matrix
        self.DVH_Data = self.DVH_Data.dropna(axis = 1) # to reduce the column (patient #) from matrix
        self.DTH_Data = self.DTH_Data.dropna(axis = 1) # to reduce the column (patient #) from matrix

        ## To further extract DTH data
        DTH_new_columns = self.further_clean_up_for_DTH(wplt = wplt) # if want to output fig,please change it to True
        self.DTH_Data = self.DTH_Data[DTH_new_columns]
        
        print("DTH Data Shape:{} and DVH Data Shape:{}".format(self.DTH_Data.shape,self.DVH_Data.shape))
        
        # save it to folder
        print("Now save clean up DVH and DTH to fixed folder!")
        self.DTH_Data.to_csv(os.path.join(self.DVH_DTH_data_path,self.strt_name+"_DTH_FurtherCleanedUp_ForSDAE.csv"))
        self.DVH_Data.to_csv(os.path.join(self.DVH_DTH_data_path,self.strt_name+"_DVH_FurtherCleanedUp_ForSDAE.csv"))

        ## standardization of DVH and DTH
        # DVH data
        self.SC_DVH = StandardScaler()
        # Norm_DVH = self.SC_DVH.fit(self.DVH_Data.T)
        Norm_DVH = self.SC_DVH.fit_transform(self.DVH_Data.T)
        self.PCA_DVH = PCA(n_components=num_components)
        self.principalComponents_DVH = self.PCA_DVH.fit_transform(Norm_DVH)
        self.principal_DVH_DATA = pd.DataFrame(data = self.principalComponents_DVH,
                                          columns = ["PC_DVH"+str(i+1) for i in range(num_components)])
        
        # DTH data
        Norm_DTH = StandardScaler().fit_transform(self.DTH_Data.T)
        PCA_DTH = PCA(n_components=num_components)
        principalComponents_DTH = PCA_DTH.fit_transform(Norm_DTH)
        self.principal_DTH_DATA = pd.DataFrame(data = principalComponents_DTH,
                                          columns = ["PC_DTH"+str(i+1) for i in range(num_components)])

        # mark the index of the case
#        principal_DVH_DATA.index = pt_id
#        principal_DTH_DATA.index = pt_id
        self.principal_DVH_DATA.index = self.DVH_Data.columns
        
#        self.principal_DTH_DATA.index = self.DTH_Data.columns
        self.principal_DTH_DATA.index = [item.split('__')[0] for item in self.DTH_Data.columns]
        
#        return self.principal_DTH_DATA, self.principal_DVH_DATA

    def Training_VS_Validation_Set(self,rtss_Volume_path,structure_name):
        '''
           shuffle the training data set and testing data set
           and determine how many percent is training set and how many 
           percent is validation set
           
           Argv:
               1.rtss_Volume is the path store volume of each rt structure 
               2.structure_name is "Parotid_L"
           
        '''
        import random
        import pandas as pd
        
        # select volume columns
        ss_volume = pd.read_csv(rtss_Volume_path)
        indexes = ss_volume.iloc[:,0]
        SS_volume = ss_volume.iloc[:,1:]
        SS_volume.index = indexes
        
#        name_index = self.principal_DVH_DATA.index
#        SS_volume = SS_volume[name_index]
        self.Target_name = 'PTV_Combinations'
        self.SS_volume = SS_volume[[structure_name,self.Target_name]]
        # combine PC and SS_Volume
        
        
        patient_names= self.principal_DVH_DATA.index
        Validat_index_DVH = random.sample(list(patient_names),3)
        Validat_index_DTH = [item for item in self.principal_DTH_DATA.index for jtem in Validat_index_DVH if jtem in item]
        
        self.Train_DVH_DATA = self.principal_DVH_DATA.drop(Validat_index_DVH)
        self.Train_DTH_DATA = self.principal_DTH_DATA.drop(Validat_index_DTH)
        
        self.Valid_DVH_DATA = self.principal_DVH_DATA.loc[Validat_index_DVH]
        self.Valid_DTH_DATA = self.principal_DTH_DATA.loc[Validat_index_DTH]
        
        self.Train_Volume = self.SS_volume.drop(Validat_index_DVH)
        self.Valid_Volume = self.SS_volume.loc[Validat_index_DVH]
        
        # combine it with volume
        self.Train_DTH_DATA_ = pd.concat([self.Train_DTH_DATA,self.Train_Volume],axis=1,sort=True)
        self.Valid_DTH_DATA_ = pd.concat([self.Valid_DTH_DATA,self.Valid_Volume],axis=1,sort=True)
        
        
#        return self.Train_DVH_DATA,self.Train_DTH_DATA_,self.Valid_DVH_DATA,self.Valid_DTH_DATA_

    def SVR_Correlation_Single_output_version(self,epsilon_1,epsilon_2,epsilon_3,epsilon_4,epsilon_5):

        '''
           This script was used for support vector regression to predict
           new DVH principle components for further analysis.
           
        '''
        from sklearn.preprocessing import StandardScaler
        from sklearn.svm import SVR
        import pandas as pd
        import numpy as np
        
        sc_X_train = StandardScaler()
        X = sc_X_train.fit_transform(self.Train_DTH_DATA_)
        sc_X_valid = StandardScaler()
        X_valid = sc_X_valid.fit_transform(self.Valid_DTH_DATA_) 
        
        sc_Y_train = StandardScaler()
        Y = sc_Y_train.fit_transform(self.Train_DVH_DATA)
        self.sc_y = StandardScaler()
        Y_valid = self.sc_y.fit_transform(self.Valid_DVH_DATA) 
        
        
        
        ######### Train model Part ###########################################
#        # Principal Component 1 #
#        PC1 = pd.DataFrame(self.Train_DVH_DATA['PC_DVH1'])
#        Y1 = sc_y1.fit_transform(PC1) # standardization
        abs_mean_error = {'PC1':[],'PC2':[],'PC3':[],'PC4':[],'PC5':[]}
        
        #PC1
        for i,item in enumerate(epsilon_1):
            Y1 = Y[:,0]  ## PC1
            PC1_regressor = SVR(kernel='linear',epsilon=item)
            PC1_regressor.fit(X,Y1)
            
            Y_predicted_PC1 = PC1_regressor.predict(X_valid)
            abs_mean_error['PC1'].append(np.mean(np.abs(Y_predicted_PC1-Y_valid[:,0])))
        
        value = min(abs_mean_error['PC1'])
        index_1 = abs_mean_error['PC1'].index(value)
        print("PC1 min absovalidation error:{} and epsilon1:{}\n".format(value,epsilon_1[index_1]))
        
        #PC2
        for i,item in enumerate(epsilon_2):
            Y2 = Y[:,1]  ## PC2
            PC2_regressor = SVR(kernel='linear',epsilon=item)
            PC2_regressor.fit(X,Y2)
            
            Y_predicted_PC2 = PC2_regressor.predict(X_valid)
            abs_mean_error['PC2'].append(np.mean(np.abs(Y_predicted_PC2-Y_valid[:,1])))
        
        value = min(abs_mean_error['PC2'])
        index_2 = abs_mean_error['PC2'].index(value)
        print("PC2 min absovalidation error:{} and epsilon2:{}\n".format(value,epsilon_2[index_2]))       
        
        #PC3
        for i,item in enumerate(epsilon_3):
            Y3 = Y[:,2]  ## PC3
            PC3_regressor = SVR(kernel='linear',epsilon=item)
            PC3_regressor.fit(X,Y3)
            
            Y_predicted_PC3 = PC3_regressor.predict(X_valid)
            abs_mean_error['PC3'].append(np.mean(np.abs(Y_predicted_PC3-Y_valid[:,2])))
        
        value = min(abs_mean_error['PC3'])
        index_3 = abs_mean_error['PC3'].index(value)
        print("PC3 min absovalidation error:{} and epsilon3:{}\n".format(value,epsilon_3[index_3]))             
        
        #PC4
        for i,item in enumerate(epsilon_4):
            Y4 = Y[:,3]  ## PC4
            PC4_regressor = SVR(kernel='linear',epsilon=item)
            PC4_regressor.fit(X,Y4)
            
            Y_predicted_PC4 = PC4_regressor.predict(X_valid)
            abs_mean_error['PC4'].append(np.mean(np.abs(Y_predicted_PC4-Y_valid[:,3])))
        
        value = min(abs_mean_error['PC4'])
        index_4 = abs_mean_error['PC4'].index(value)
        print("PC4 min absovalidation error:{} and epsilon4:{}\n".format(value,epsilon_4[index_4]))               
        
        #PC5
        for i,item in enumerate(epsilon_5):
            Y5 = Y[:,4]  ## PC5
            PC5_regressor = SVR(kernel='linear',epsilon=item)
            PC5_regressor.fit(X,Y5)
            
            Y_predicted_PC5 = PC5_regressor.predict(X_valid)
            abs_mean_error['PC5'].append(np.mean(np.abs(Y_predicted_PC5-Y_valid[:,4])))
        
        value = min(abs_mean_error['PC5'])
        index_5 = abs_mean_error['PC5'].index(value)
        print("PC5 min absovalidation error:{} and epsilon5:{}\n".format(value,epsilon_5[index_5]))  
        
        # Predicted PC components 
        Y_predicted_PC = []
        epsilon = [epsilon_1[index_1],epsilon_2[index_2],epsilon_3[index_3],epsilon_4[index_4],epsilon_5[index_5]]
        for kk,epsilon_ in enumerate(epsilon):
            
            PC_regressor = SVR(kernel='linear',epsilon=epsilon_)
            PC_regressor.fit(X,Y[:,kk])
            Y_predicted_PC.append(PC_regressor.predict(X_valid))
        
        ## predicted Y
        self.Y_predicted = np.array(Y_predicted_PC).T
        self.y_pred = self.Y_predicted
        ## inverse transform
        self.y_test = Y_valid

        self.y_test_ = pd.DataFrame(self.y_test,index = self.Valid_DVH_DATA.index,columns = self.Valid_DVH_DATA.columns)
        self.y_pred_ = pd.DataFrame(self.y_pred,index = self.Valid_DVH_DATA.index,columns = self.Valid_DVH_DATA.columns)
        
        
        return self.y_pred_, self.y_test_


    def SVR_Correlation_Multi_output_version(self,eplison_):

        '''
           This script was used for support vector regression to predict
           new DVH principle components for further analysis.
           
           Args:
               1. eplison_ is the sensitivity factor for SVR model 
        '''
        from sklearn.preprocessing import StandardScaler
        from sklearn.datasets import make_regression
        from sklearn.multioutput import MultiOutputRegressor
        from sklearn.svm import SVR
        from sklearn.model_selection import train_test_split
        from sklearn.metrics import mean_squared_error, mean_absolute_error
        import pandas as pd
        
        sc_X = StandardScaler()
        self.sc_y = StandardScaler()
        X_train = sc_X.fit_transform(self.Train_DTH_DATA_)
        X_valid = sc_X.fit_transform(self.Valid_DTH_DATA_) 
        y_train = self.sc_y.fit_transform(self.Train_DVH_DATA)
        self.y_test = self.sc_y.fit_transform(self.Valid_DVH_DATA)
        ######### Train model Part ###########################################
        # Create the SVR regressor
        svr = SVR(epsilon=eplison_)
        
        # Create the Multioutput Regressor
        mor = MultiOutputRegressor(svr)
        
        # Train the regressor
        mor = mor.fit(X_train, y_train)
        
        # Generate predictions for testing data
        self.y_pred = mor.predict(X_valid)
        
        # Evaluate the regressor
        mse_1 = mean_squared_error(self.y_test[:,0], self.y_pred[:,0])
        mse_2 = mean_squared_error(self.y_test[:,1], self.y_pred[:,1])
        mse_3 = mean_squared_error(self.y_test[:,2], self.y_pred[:,2])
        mse_4 = mean_squared_error(self.y_test[:,3], self.y_pred[:,3])
        mse_5 = mean_squared_error(self.y_test[:,4], self.y_pred[:,4])
        print('MSE for first regressor: {}'.format(mse_1))
        print('MSE for second regressor: {}'.format(mse_2))
        print('MSE for third regressor: {}'.format(mse_3))
        print('MSE for forth regressor: {}'.format(mse_4))
        print('MSE for fifth regressor: {}'.format(mse_5))

        mae_1 = mean_absolute_error(self.y_test[:,0], self.y_pred[:,0])
        mae_2 = mean_absolute_error(self.y_test[:,1], self.y_pred[:,1])
        mae_3 = mean_absolute_error(self.y_test[:,2], self.y_pred[:,2])
        mae_4 = mean_absolute_error(self.y_test[:,3], self.y_pred[:,3])
        mae_5 = mean_absolute_error(self.y_test[:,4], self.y_pred[:,4])
        print('MAE for first regressor: {}'.format(mae_1))
        print('MAE for second regressor: {}'.format(mae_2))
        print('MAE for third regressor: {}'.format(mae_3))
        print('MAE for forth regressor: {}'.format(mae_4))
        print('MAE for fifth regressor: {}'.format(mae_5))
        
        self.y_test_ = pd.DataFrame(self.y_test,index = self.Valid_DVH_DATA.index,columns = self.Valid_DVH_DATA.columns)
        self.y_pred_ = pd.DataFrame(self.y_pred,index = self.Valid_DVH_DATA.index,columns = self.Valid_DVH_DATA.columns)
        
        return self.y_test_, self.y_pred_


    def DVH_Reconstruct(self):
        '''
           This script was used for reconstruct the DVH from the predicted 
           DVH principle components.
        '''
#        from sklearn.preprocessing import StandardScaler
        import numpy as np
        import pandas as pd
        # to list the validated DVH 
        pt_names = list(self.Valid_DVH_DATA.index)
        DVH_Valid_Truth = self.DVH_Data[pt_names]
        
        # reconstruct the predicted DVH
        components = self.sc_y.inverse_transform(self.y_pred)
        self.DVH_Valid_Pred = self.SC_DVH.inverse_transform(np.dot(components,self.PCA_DVH.components_))
        self.DVH_Valid_Pred[self.DVH_Valid_Pred > 100] = 100 # to cut out the high volume percent
        self.DVH_Valid_Pred[self.DVH_Valid_Pred < 0] = 0 # to cut out the volume value below zero

        new_columns = [item+"_PCA_SVR_Pred" for item in DVH_Valid_Truth.columns]
        self.DVH_Valid_Pred = pd.DataFrame(self.DVH_Valid_Pred.T,index =DVH_Valid_Truth.index, columns = new_columns)
        
        return self.DVH_Valid_Pred, DVH_Valid_Truth

class SDAE_FCNN_Model:
    '''
       This class was followed by a research paper to use SADE(Stack_Encode_Decode)
       framework to reduce the data dimension of DTH and DVH which is similar to PCA 
       methods and Fully Connected Neural Network was used for correlate the DTH with DVH for further 
       prediction of DVH.

       The data preprocessing for SDAE is similar to PCA+SVR while it should be noted, SDAE model should be carefully
       constructed to make it perform better and efficient.

    '''
    def __init__(self,DVH_DTH_data_path,strt_name):
        '''
          Args:
          1. DVH_DTH_data_path is a place for storing DTH and DVH data.
          2. strt_name is a variable for which structure need to be dealt with.
        '''
        self.DVH_DTH_data_path = DVH_DTH_data_path
        self.strt_name = strt_name

    def data_clean_preparation(self,compress,file_name_DVH,file_name_DTH):
        '''
          This script was used for divide the data into train and valid set
          Args:
          compress: to reduce the original data index spacing of DTH and DVH
          file_name_DVH: "DVH_FurtherCleanedUp_ForSDAE.csv"
          file_name_DTH: "DTH_FurtherCleanedUp_ForSDAE.csv"
        '''
        import pandas as pd
#        import numpy as np
        import os
        import random
        from sklearn.preprocessing import StandardScaler

        # read data into RAM
        self.DVH_Data = pd.read_csv(os.path.join(self.DVH_DTH_data_path,file_name_DVH))
        self.DTH_Data = pd.read_csv(os.path.join(self.DVH_DTH_data_path,file_name_DTH))

        DVH_index = self.DVH_Data.iloc[:,0]
        DTH_index = self.DTH_Data.iloc[:,0]

        self.DVH_Data = self.DVH_Data.iloc[:,1:]
        self.DTH_Data = self.DTH_Data.iloc[:,1:]

        self.DVH_Data.index = DVH_index
        self.DTH_Data.index = DTH_index

        self.Validat_index_DVH = random.sample(list(self.DVH_Data.columns),5)
        
#        self.Validat_index_DVH = random.sample(list(self.DVH_Data.columns),5)
        # Validat_index_DVH = ['HNSCC-01-0161_rtss','HNSCC-01-0075_rtss','HNSCC-01-0053_rtss','HNSCC-01-0045']
#        self.Validat_index_DTH = [item for item in self.DTH_Data.columns for jtem in self.Validat_index_DVH if jtem in item]

#        self.Valid_DVH_DATA = self.DVH_Data[self.Validat_index_DVH]
#        self.Valid_DTH_DATA = self.DTH_Data[self.Validat_index_DTH]
#
#        self.Train_DVH_DATA = self.DVH_Data.drop(self.Validat_index_DVH,axis=1)
#        self.Train_DTH_DATA = self.DTH_Data.drop(self.Validat_index_DTH,axis=1)
        
        ## compress dimensions from real to 100 and do data standardization
        self.SC1 = StandardScaler()
        self.DVH_Data_ = self.SC1.fit_transform(self.DVH_Data.iloc[::compress,:].T)
        self.SC2 = StandardScaler()
        self.DTH_Data_ = self.SC2.fit_transform(self.DTH_Data.iloc[::compress,:].T)
#        self.SC3 = StandardScaler()
#        self.Train_DVH_DATA_ = self.SC3.fit_transform(self.Train_DVH_DATA.iloc[::compress,:].T)
#        self.SC4 = StandardScaler()
#        self.Train_DTH_DATA_ = self.SC4.fit_transform(self.Train_DTH_DATA.iloc[::compress,:].T)
#        self.SC5 = StandardScaler()
#        self.Valid_DVH_DATA_ = self.SC5.fit_transform(self.Valid_DVH_DATA.iloc[::compress,:].T)
#        self.SC6 = StandardScaler()
#        self.Valid_DTH_DATA_ = self.SC6.fit_transform(self.Valid_DTH_DATA.iloc[::compress,:].T)
        self.SC3 = StandardScaler()
        test_data = self.DVH_Data[self.Validat_index_DVH]
        self.DVH_DATA_train = self.SC3.fit_transform(test_data.iloc[::compress,:].T)

        return self.SC1, self.SC2, self.SC3

    def SDAE_Compression(self,type,epoch):
        '''
          This script was used for data compression from high dimensions
          to low dimensios(5) with the Encoder-To-Decoder

          Args:
          1. type has two selection: DVH and DTH for different data compress
        '''
        import tensorflow.keras as keras
        from tensorflow.keras import layers
        import pandas as pd
#        from tensorflow.keras.callbacks import TensorBoard

        if type == "DVH":
            input_data_type = self.DVH_Data_
            x_train = self.DVH_Data_
#            x_test = self.Valid_DVH_DATA_
        else:
            input_data_type = self.DTH_Data_
            x_train = self.DTH_Data_
#            x_test = self.Valid_DTH_DATA_

        # This is the size of our encoded representations
        encoding_dim = 10 # to compress the data to 5 dimensions

        # This is our input image
        input_data = keras.Input(shape=(input_data_type.shape[1],)) # for DVH 8142, for DTH 18229

        encoded = layers.Dense(64, activation='relu')(input_data)
        encoded = layers.Dense(128, activation='relu')(encoded)
        encoded = layers.Dense(256, activation='relu')(encoded)
        encoded = layers.Dense(128, activation='relu')(encoded)
        encoded = layers.Dense(64, activation='relu')(encoded)
        encoded = layers.Dense(32, activation='relu')(encoded)
        encoded = layers.Dense(16, activation='relu')(encoded)
        encoded = layers.Dense(encoding_dim, activation='relu')(encoded)

        decoded = layers.Dense(16, activation='relu')(encoded)
        decoded = layers.Dense(32, activation='relu')(decoded)
        decoded = layers.Dense(64, activation='relu')(decoded)
        decoded = layers.Dense(128, activation='relu')(decoded)
        decoded = layers.Dense(256, activation='relu')(decoded)
        decoded = layers.Dense(128, activation='relu')(decoded)
        decoded = layers.Dense(64, activation='relu')(decoded)
        decoded = layers.Dense(input_data_type.shape[1], activation='sigmoid')(decoded)


        self.autoencoder = keras.Model(input_data, decoded)
        # mse = keras.losses.MeanSquaredError()
        self.autoencoder.compile(optimizer='adam', loss='mean_squared_error')

        
        self.autoencoder.fit(x_train, x_train,
                epochs=epoch,
                batch_size=128,
                shuffle=True,
                validation_split = 0.3
                )
        
        
        self.encoder = keras.Model(input_data, encoded) # for output the compressed data

        # This is our encoded (10-dimensional) input
        self.encoded_input = keras.Input(shape=(encoding_dim,))
        # Retrieve the last layer of the autoencoder model
        deco = self.autoencoder.layers[-8](self.encoded_input)
        deco = self.autoencoder.layers[-7](deco)
        deco = self.autoencoder.layers[-6](deco)
        deco = self.autoencoder.layers[-5](deco)
        deco = self.autoencoder.layers[-4](deco)
        deco = self.autoencoder.layers[-3](deco)
        deco = self.autoencoder.layers[-2](deco)
        self.deco = self.autoencoder.layers[-1](deco)
        # Create the decoder model
        self.decoder = keras.Model(self.encoded_input, self.deco)


        self.compressed_train_data = self.encoder.predict(x_train)  ## dimension 5
#        self.compressed_test_data = self.encoder.predict(x_test)
        self.decoder_train_data = self.autoencoder.predict(x_train)  ## original dimensions
#        self.decoder_test_data = self.autoencoder.predict(x_test)

        self.compressed_train_data = pd.DataFrame(self.compressed_train_data,index = self.DVH_Data.columns)
        self.decoder_train_data    = pd.DataFrame(self.decoder_train_data,index = self.DVH_Data.columns)

        return self.compressed_train_data.T,self.decoder_train_data.T

    def _Data_preparation_1DCNN(self,compress,epochs):
        '''
           This function was used for data preparation for 1DCNN
        '''
        import os
#        import random 
        import pandas as pd


#        self.Validat_index_ = random.sample(list(self.DVH_Data.columns),5)
        # compress = 80
        strt_name = ["Parotid_L","Parotid_R","SpinalCord","BrainStem"]
        self.SC_Store,self.DATA_Prepared = {},{}
        self.decoder_machine = {}
        ## Prepare the training data (ParotidL4,ParotidR4,BrainStem4,SpinalCord4,Vol_ParotidL,Vol_ParotidR,Vol_BS,Vol_SC)
        ## Input data should be matrix has dimensions: (patient_num, 20)
        ## output data should be matrix has dimensions: (patient_num, 5)
#        self.Validat_index_DVH = random.sample(list(self.DVH_Data.columns),5)
        self.Validat_index_DTH = [item for item in self.DTH_Data.columns for jtem in self.Validat_index_DVH if jtem in item]
        for name in strt_name:
            S1,S2,S3 = self.data_clean_preparation(compress,
                                                name+"_DVH_FurtherCleanedUp_ForSDAE.csv",
                                                name+"_DTH_FurtherCleanedUp_ForSDAE.csv")
#            self.SC_Store[name] = {"S1":S1,"S2":S2,"S3":S3}
            if name == self.strt_name:
                
                compressed_data_DTH,decoder_data_DTH = self.SDAE_Compression("DTH",epochs)   # to compress data firstly
                self.DATA_Prepared[name+"_DTH_test"] = compressed_data_DTH[self.Validat_index_DVH]
                self.DATA_Prepared[name+"_DTH_train"] = compressed_data_DTH.drop(self.Validat_index_DVH,axis=1)
                
                compressed_data_DVH,decoder_data_DVH = self.SDAE_Compression("DVH",epochs)
                self.DATA_Prepared[name+"_DVH_test"] = compressed_data_DVH[self.Validat_index_DVH]
                self.DATA_Prepared[name+"_DVH_train"] = compressed_data_DVH.drop(self.Validat_index_DVH,axis=1)
                self.decoder_machine[name] = self.decoder
                
            else:
                
                compressed_data_DTH,decoder_data_DTH = self.SDAE_Compression("DTH",epochs)   # to compress data firstly               
    
                self.DATA_Prepared[name+"_DTH_test"] = compressed_data_DTH[self.Validat_index_DVH]
                self.DATA_Prepared[name+"_DTH_train"] = compressed_data_DTH.drop(self.Validat_index_DVH,axis=1)

        Volume = pd.read_csv(os.path.join(self.DVH_DTH_data_path,"Volume_Clean_New.csv"))
        index_name = Volume.iloc[:,0]
        Volume = Volume.iloc[:,1:]
        Volume.index = index_name
        Volume = Volume.T

        drop_index = self.Validat_index_DVH
        # Validat_index_DVH = random.sample(list(self.DVH_Data.columns),10)
        self.test_Volume = Volume[drop_index]
        self.train_Volume = Volume.drop(columns=drop_index)
        
        
        ### 
        
        

    def OneD_CNN_Correlation(self,Model_Type_For_DVH,CNN1D_epoch):
        '''
          This work was used for building up a 1D-CNN to help correlate 
          DTH and DVH (Parotid L, Parotid R, Spinal Cord, Brain Stem) DTH and DVH

          Model_Type_For_DVH : "Parotid_L,Parotid_R,SpinalCord,BrainStem"
        '''
        from tensorflow.keras.models import Sequential
        from tensorflow.keras.layers import Dense, Conv1D, Flatten, MaxPooling1D, Dropout
        from tensorflow.keras.optimizers import RMSprop
        from tensorflow.keras.callbacks import TensorBoard
        import pandas as pd 
        import numpy as np
        import os
        from sklearn.preprocessing import StandardScaler
        
        input_train_data = np.column_stack((self.DATA_Prepared['BrainStem_DTH_train'].T.values,
                                           self.DATA_Prepared['SpinalCord_DTH_train'].T.values,
                                           self.DATA_Prepared['Parotid_L_DTH_train'].T.values,
                                           self.DATA_Prepared['Parotid_R_DTH_train'].T.values,
                                           self.train_Volume.T.values))

        input_test_data = np.column_stack((self.DATA_Prepared['BrainStem_DTH_test'].T.values,
                                           self.DATA_Prepared['SpinalCord_DTH_test'].T.values,
                                           self.DATA_Prepared['Parotid_L_DTH_test'].T.values,
                                           self.DATA_Prepared['Parotid_R_DTH_test'].T.values,
                                           self.test_Volume.T.values))

        print("input_train_data's shape:{}, input_test_data's shape:{}".format(input_train_data.shape,input_test_data.shape))


        if Model_Type_For_DVH == 'Parotid_L' or Model_Type_For_DVH == 'Parotid_R' or Model_Type_For_DVH == 'BrainStem' or Model_Type_For_DVH == 'SpinalCord':
            y_train = self.DATA_Prepared[Model_Type_For_DVH+'_DVH_train'].T.values
            self.y_test = self.DATA_Prepared[Model_Type_For_DVH+'_DVH_test'].T.values


        ## Standardization of training data
        self.train_SC = StandardScaler()
        x_train = self.train_SC.fit_transform(input_train_data)
        self.test_SC = StandardScaler()
        x_test = self.test_SC.fit_transform(input_test_data)
        x_train = x_train.reshape((x_train.shape[0],x_train.shape[1],1))
        x_test = x_test.reshape((x_test.shape[0],x_test.shape[1],1))

        self.train_y_SC = StandardScaler()
        y_train = self.train_y_SC.fit_transform(y_train)

        self.test_y_SC = StandardScaler()
        self.y_test = self.test_y_SC.fit_transform(self.y_test)

        ## Model Build-Up 
        model = Sequential()
        model.add(Dropout(0.2, input_shape=(input_train_data.shape[1], 1)))
        model.add(Conv1D(20, kernel_size = 3, activation="relu"))
        model.add(Conv1D(32, kernel_size = 3, activation="relu"))
        model.add(MaxPooling1D(2))
        # model.add(MaxPooling1D(pool_size=2,strides=1,padding='valid'))

        model.add(Conv1D(64, kernel_size = 3, activation="relu"))
        model.add(Conv1D(128, kernel_size = 3, activation="relu"))
        model.add(MaxPooling1D(2))
        # model.add(MaxPooling1D(pool_size=2,strides=1,padding='valid'))
        model.add(Conv1D(256, kernel_size = 3, activation="relu"))
        model.add(Conv1D(512, kernel_size = 3, activation="relu"))       
        model.add(Flatten())

        model.add(Dense(128, activation="relu"))
        model.add(Dense(32, activation="relu"))
        model.add(Dense(10))
        
        model.compile(loss="mse", optimizer="adam", metrics=['accuracy'])
        print(model.summary())

        self.history = model.fit(x_train, y_train,
                  epochs = CNN1D_epoch,
                  shuffle=True,
                  batch_size=20,
                  validation_split = 0.3)
        
        self.predict_y = model.predict(x_test)

        score = model.evaluate(x_test, self.y_test, verbose=0)
        print('Test loss:', score[0])
        print('Test accuracy:', score[1])

        return self.predict_y,self.y_test


    def data_clean_preparation_inverse(self,compress,file_name_DVH,file_name_DTH,decoder_predict_y):

        '''
          This script was used for divide the data into train and valid set
          Args:
          compress: to reduce the original data index spacing of DTH and DVH
          file_name_DVH: "DVH_FurtherCleanedUp_ForSDAE.csv"
          file_name_DTH: "DTH_FurtherCleanedUp_ForSDAE.csv"
        '''
        import pandas as pd
#        import numpy as np
        import os
#        import random
        from sklearn.preprocessing import StandardScaler

        # read data into RAM
        self.DVH_Data = pd.read_csv(os.path.join(self.DVH_DTH_data_path,file_name_DVH))
        self.DTH_Data = pd.read_csv(os.path.join(self.DVH_DTH_data_path,file_name_DTH))

        DVH_index = self.DVH_Data.iloc[:,0]
        DTH_index = self.DTH_Data.iloc[:,0]

        self.DVH_Data = self.DVH_Data.iloc[:,1:]
        self.DTH_Data = self.DTH_Data.iloc[:,1:]

        self.DVH_Data.index = DVH_index
        self.DTH_Data.index = DTH_index

#        self.Validat_index_DVH = random.sample(list(self.DVH_Data.columns),5)
        # Validat_index_DVH = ['HNSCC-01-0161_rtss','HNSCC-01-0075_rtss','HNSCC-01-0053_rtss','HNSCC-01-0045']
#        self.Validat_index_DTH = [item for item in self.DTH_Data.columns for jtem in self.Validat_index_DVH if jtem in item]

#        self.Valid_DVH_DATA = self.DVH_Data[self.Validat_index_DVH]
#        self.Valid_DTH_DATA = self.DTH_Data[self.Validat_index_DTH]
#
#        self.Train_DVH_DATA = self.DVH_Data.drop(self.Validat_index_DVH,axis=1)
#        self.Train_DTH_DATA = self.DTH_Data.drop(self.Validat_index_DTH,axis=1)
        
        ## compress dimensions from real to 100 and do data standardization
#        self.SC1 = StandardScaler()
#        self.DVH_Data_ = self.SC1.fit_transform(self.DVH_Data.iloc[::compress,:].T)
#        self.SC2 = StandardScaler()
#        self.DTH_Data_ = self.SC2.fit_transform(self.DTH_Data.iloc[::compress,:].T)
#        self.SC3 = StandardScaler()
#        self.Train_DVH_DATA_ = self.SC3.fit_transform(self.Train_DVH_DATA.iloc[::compress,:].T)
#        self.SC4 = StandardScaler()
#        self.Train_DTH_DATA_ = self.SC4.fit_transform(self.Train_DTH_DATA.iloc[::compress,:].T)
#        self.SC5 = StandardScaler()
#        self.Valid_DVH_DATA_ = self.SC5.fit_transform(self.Valid_DVH_DATA.iloc[::compress,:].T)
#        self.SC6 = StandardScaler()
#        self.Valid_DTH_DATA_ = self.SC6.fit_transform(self.Valid_DTH_DATA.iloc[::compress,:].T)
        SC3 = StandardScaler()
        test_data = self.DVH_Data[self.Validat_index_DVH]
        self.DVH_DATA_train = SC3.fit_transform(test_data.iloc[::compress,:].T)
        decoder_predict_y = SC3.inverse_transform(decoder_predict_y)
        
                
        decoder_predict_y = pd.DataFrame(decoder_predict_y, index = self.Validat_index_DVH,columns=test_data.iloc[::compress,:].index)
    
        return decoder_predict_y

        

    def DVH_Reconstruct(self,strt_name):
        '''
           This script is intended to reconstruct DVH curves from predicted
           results
           
           e.g. for self.predict_y, it should firstly, inverse_standardization
           and then decoder to DVH data
           e.g. for self.y_test, it should firstly, inverse_standardization
           and then also decoder to DVH data.
        '''
#        import pandas as pd
#        from sklearn.preprocessing import StandardScaler
        compress = 100
        
        ## inverse standardization of data
        predict_y_ = self.test_y_SC.inverse_transform(self.predict_y)
#        test_y_    = self.test_y_SC.inverse_transform(self.y_test)
        
        ## decode the data for DVH evaluation
        decoder_predict_y = self.decoder_machine[strt_name].predict(predict_y_)
#        decoder_test_y    = self.decoder_machine[strt_name].predict(test_y_)

                
        ## inverse standardization
        decoder_predict_y_ = self.data_clean_preparation_inverse(compress,strt_name+"_DVH_FurtherCleanedUp_ForSDAE.csv",
                                                strt_name+"_DTH_FurtherCleanedUp_ForSDAE.csv",decoder_predict_y)
#        decoder_test_y_ = self.data_clean_preparation_inverse(compress,strt_name+"_DVH_FurtherCleanedUp_ForSDAE.csv",
#                                                strt_name+"_DTH_FurtherCleanedUp_ForSDAE.csv",decoder_test_y)


        
        return decoder_predict_y_

