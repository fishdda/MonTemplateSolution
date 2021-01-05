class DATA_CLEAN:

    def __init__(self, DICOM_name, RTP_name, goal_path, TPS_name):
        self.DICOM_name = DICOM_name
        self.RTP_name = RTP_name
        self.goal_path = goal_path
        self.TPS_name = TPS_name

    def extract_DICOM(self):
        '''
           This function was used for extracting DICOM file
        '''

        import pydicom 
        import os
        dcm = pydicom.read_file(self.DICOM_name,force=True)
        # TPS = 'Monaco'
        mu_density = {'Gantry':[],'MLC':[],'JAW':[],'MU':[],'Absolute_MU':[]}

        for i in range(len(dcm.BeamSequence)):
            
            mu_density['Absolute_MU'].append(float(dcm.FractionGroupSequence[0].ReferencedBeamSequence[i].BeamMeterset))
            
            for j in range(len(dcm.BeamSequence[i].ControlPointSequence)):
                
                mu_density['Gantry'].append(int(dcm.BeamSequence[i].ControlPointSequence[0].GantryAngle))

                if self.TPS_name == 'Pinnacle':
        #             print(i,j)
                    MU = float(dcm.BeamSequence[i].ControlPointSequence[j].CumulativeMetersetWeight)*float(mu_density['Absolute_MU'][i])/100
        #             print('MU is {}'.format(MU))
                    mu_density['MU'].append(MU)
                    if len(dcm.BeamSequence[i].ControlPointSequence[j].BeamLimitingDevicePositionSequence) == 3:
                        mu_density['MLC'].append(dcm.BeamSequence[i].ControlPointSequence[j].BeamLimitingDevicePositionSequence[2].LeafJawPositions)
                        mu_density['JAW'].append(dcm.BeamSequence[i].ControlPointSequence[j].BeamLimitingDevicePositionSequence[1].LeafJawPositions)
                    else:
                        mu_density['MLC'].append(dcm.BeamSequence[i].ControlPointSequence[j].BeamLimitingDevicePositionSequence[0].LeafJawPositions)
                        mu_density['JAW'].append([None,None])
                        print('it has no jaw positions')
                        print('This is {}th beam {}th control point'.format(i+1,j+1))
                        print('MU is {}'.format(MU))
                elif self.TPS_name == 'Monaco':
                    MU = float(dcm.BeamSequence[i].ControlPointSequence[j].CumulativeMetersetWeight)
                    mu_density['MU'].append(MU)
                    mu_density['MLC'].append(dcm.BeamSequence[i].ControlPointSequence[j].BeamLimitingDevicePositionSequence[1].LeafJawPositions)
                    mu_density['JAW'].append(dcm.BeamSequence[i].ControlPointSequence[j].BeamLimitingDevicePositionSequence[0].LeafJawPositions)
                    

        #further clearning
        PLAN = {key:[] for key in mu_density.keys()}
        PLAN['MU_absolute'] = mu_density['Absolute_MU']
        for i,item in enumerate(mu_density['MU']):
            if i==0:
                PLAN['JAW'].append(mu_density['JAW'][i])
                PLAN['MU'].append(mu_density['MU'][i])
                PLAN['Gantry'].append(mu_density['Gantry'][i])
                PLAN['MLC'].append(mu_density['MLC'][i])
            else:
                if mu_density['MU'][i] != mu_density['MU'][i-1]:
                    PLAN['JAW'].append(mu_density['JAW'][i])
                    PLAN['MU'].append(mu_density['MU'][i])
                    PLAN['Gantry'].append(mu_density['Gantry'][i])
                    PLAN['MLC'].append(mu_density['MLC'][i])

        self.PLAN = PLAN
        return self.PLAN   

    def extract_RTP(self):
        '''
           This function was used for extracting RTP file exported from MOSAIQ
        '''
        line,pointer = [],[]
        absolute_MU = []
        count = -1
        RTP = {'Gantry':[],'MLC':[],'JAW':[],'MU':[],'MU_absolute':[]}
        with open(self.RTP_name, "r+", encoding = "ISO-8859-1") as f:
            line1 = f.readline()    
            line.append(line1)
            while line1:        
                pointer.append(f.tell())  #record the pointer loaction to help write        
                line1 = f.readline()        
                line.append(line1)
        for i,item in enumerate(line[2:-1]):
            TempLine = item.split('"')  
            if TempLine[1] == 'CONTROL_PT_DEF': # to ensure the head
                RTP['Gantry'].append(float(TempLine[27].split(' ')[-1]))
                PP_ = []
                for item1 in TempLine:
                    if ' ' in item1:
                        PP_.append(item1.split(' ')[-1])
                    else:
                        PP_.append(item1)

        #         print('This is {}th control point'.format(i))
                PP_ = [item for item in PP_ if item != ',' and item != '']
                if self.TPS_name == 'Pinnacle':
                    MLC = PP_[27:-2]
        #             print(len(MLC))
                    MLC_A = [float(item) for item in MLC[0:80]]
                    MLC_B = [float(item) for item in MLC[80:]]
                    RTP['MLC'].append(MLC_A+MLC_B)
                    RTP['MU'].append(float(PP_[7])*RTP['MU_absolute'][count])
                    JAW1 = float(TempLine[47].split(' ')[-1])
                    JAW2 = float(TempLine[49].split(' ')[-1])
                    RTP['JAW'].append([JAW1,JAW2])
                elif self.TPS_name == 'Monaco':
                    MLC = PP_[30:-2]
        #             print(len(MLC))
                    MLC_A = [float(item) for item in MLC[0:80]]
                    MLC_B = [float(item) for item in MLC[80:]]
                    RTP['MLC'].append(MLC_A+MLC_B)
                    RTP['MU'].append(float(PP_[7]))
                    JAW1 = float(TempLine[47].split(' ')[-1])
                    JAW2 = float(TempLine[49].split(' ')[-1])
                    RTP['JAW'].append([JAW1,JAW2])
            elif TempLine[1] == "FIELD_DEF":
                RTP['MU_absolute'].append(float(TempLine[13].split(' ')[-1]))
                count = count + 1

        # further cleaning for RTP information
        RTP_ = {key:[] for key in RTP.keys()}
        RTP_['MU_absolute'] = RTP['MU_absolute']
        for i,item in enumerate(RTP['MU']):
            if i==0:
                RTP_['JAW'].append(RTP['JAW'][i])
                RTP_['MU'].append(RTP['MU'][i])
                RTP_['Gantry'].append(RTP['Gantry'][i])
                RTP_['MLC'].append(RTP['MLC'][i])
            else:
                if RTP['MU'][i] != RTP['MU'][i-1]:
                    RTP_['JAW'].append(RTP['JAW'][i])
                    RTP_['MU'].append(RTP['MU'][i])
                    RTP_['Gantry'].append(RTP['Gantry'][i])
                    RTP_['MLC'].append(RTP['MLC'][i])

        self.RTP = RTP_
        return self.RTP
        
    def CHECKING(self):
        '''
          Check all the control points including gantry angles, plan MU, MU difference, MLC
          difference, JAW difference etc.

          Two tables would be exported, summary and MLC details
        '''

        import numpy as np
        import pandas as pd
        import os
        tolerance_mlc = 0.1 # 0.1mm MLC position tolerance 
        tolerance_gantry = 0.1 # degree
        tolerance_jaw = 0.1 # 0.1mm
        tolerance_mu = 0.1 # 

        MLC_err = np.zeros([len(self.RTP['Gantry']),160])
        for i in range(len(self.PLAN['Gantry'])):
            MLC_err[i,:] = np.round(np.array(self.PLAN['MLC'][i]) - np.array(self.RTP['MLC'][i])*10)
        columns = ['Leaf A #'+ str(i+1) + ' position error(mm)' for i in range(80)] + ['Leaf B #'+ str(i+1) + ' position error(mm)' for i in range(80)]
        MLC_err_table = pd.DataFrame(MLC_err,columns = columns)
        MLC_state = []
        for i in range(len(self.PLAN['Gantry'])):
            X = MLC_err_table.iloc[i]<tolerance_mlc
            MLC_state.append(X.all())

        MLC_err_table.to_csv(os.path.join(self.goal_path,'MLC errors deatils.csv'), index = True)  # output
        print("successfully export MLC errors details!")

        col = ['Control Point #',
            'Gantry Angle(PLAN)','Gantry Angle(RTP)',
            'Plan MU','RTP MU',
            'Plan JAW1(mm)','Plan JAW2(mm)','RTP JAW1(mm)','RTP JAW2(mm)',
            'Gantry Angle Difference(1.0 means pass)',
            'MU difference (1.0 means pass)',
            'JAW difference(1.0 means pass)',
            'MLC difference(1.0 means pass, 0 means errors exist)']
        data = np.zeros([len(self.RTP['Gantry']),13])
        for i in range(len(self.PLAN['Gantry'])): 
            data[i,0] = int(i+1)
            S = np.array(self.PLAN['MLC'][i]) == np.array(self.RTP['MLC'][i])*10
            S2 = np.array(self.PLAN['JAW'][i]) == np.array(self.RTP['JAW'][i])*10
            # Gantry difference checking
            data[i,1] = self.PLAN['Gantry'][i]
            data[i,2] = self.RTP['Gantry'][i]
        #     data[i,3] = np.round(self.PLAN['Gantry'][i] - RTP['Gantry'][i]) 
            data[i,9] = (self.PLAN['Gantry'][i] - self.RTP['Gantry'][i]) < tolerance_gantry
            # MU difference checking 
            data[i,3] = np.array(self.PLAN['MU'][i])
            data[i,4] = np.array(self.RTP['MU'][i])
        #     data[i,7] = np.round(np.array(self.PLAN['MU'][i])-np.array(RTP['MU'][i]))
            data[i,10] = (np.array(self.PLAN['MU'][i])-np.array(self.RTP['MU'][i])) < tolerance_mu
            # JAW position checking including JAW1, and JAW2
            data[i,5] = np.array(self.PLAN['JAW'][i])[0]
            data[i,6] = np.array(self.PLAN['JAW'][i])[1]
            JAW1 = np.array(self.RTP['JAW'][i])*10
            data[i,7] = JAW1[0]
            data[i,8] = JAW1[1]
            data[i,11] = S2.all()
            # MLC position checking
            data[i,12] = MLC_state[i]
            
        df = pd.DataFrame(data,columns = col)   

        df.to_csv(os.path.join(self.goal_path,'CheckingList.csv'), index = True) 
        print("successfully export CheckingList!")

