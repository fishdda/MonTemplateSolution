class DATAPROCESS:
    '''
    Currently, this class is only specific to VMAT plans 
    1. transfer all the dicom data to npy
    2. generate mu density from dicom and trf, compare them together.
    '''

    def __init__(self,RAW_DATA_PATH,NEW_DATA_PATH):
        self.RAW_DATA_PATH = RAW_DATA_PATH
        self.NEW_DATA_PATH = NEW_DATA_PATH
        self.PLAN = {}
        self.Mapping = {}
        self.DVHSTAT = {}
        self.DVHDATA = {}

    def TRANSFER_RAW_DATA(self):

        import numpy as np 
        import os
        import pydicom 
        import json


        # deal with CT image data, dose data and MLC & MU intensity

        file_names = os.listdir(self.RAW_DATA_PATH)

        for j,file_nam in enumerate(file_names):
            CT_ = []
            print('This is {}th patient'.format(j))
            temp_files = os.listdir(os.path.join(self.RAW_DATA_PATH,file_nam)) # all CT names
            for dcm_nam in temp_files:
                # 
                dcm = pydicom.read_file(os.path.join(self.RAW_DATA_PATH,file_nam,dcm_nam),force=True)
                if dcm.Modality == 'CT':
                    # print(dcm.Modality)
                    dcm.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
                    img = dcm.pixel_array

                    CT_.append(img)
                    intercept = dcm.RescaleIntercept
                    slope = dcm.RescaleSlope

                # The intercept is usually -1024, so air is approximately 0
                elif dcm.Modality == 'RTSTRUCT':
                    

                    for i in range(len(dcm.StructureSetROISequence)):
                        self.Mapping[dcm.StructureSetROISequence[i].ROINumber] = dcm.StructureSetROISequence[i].ROIName

                elif dcm.Modality == 'RTDOSE' and dcm.DoseSummationType == 'PLAN':

                    dcm.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
                    DOSE = np.array(dcm.pixel_array * dcm.DoseGridScaling)
                    print('The DOSE array shape:{}'.format(DOSE.shape))
                    np.save(self.NEW_DATA_PATH + file_nam + '_Dose.npy',DOSE)

                    for i in range(len(dcm.DVHSequence)):

                        ROI_name = self.Mapping[dcm.DVHSequence[i].DVHReferencedROISequence[0].ReferencedROINumber]
                        self.DVHSTAT[ROI_name] = {'MeanDose':dcm.DVHSequence[i].DVHMeanDose,
                                                  'MaxDose':dcm.DVHSequence[i].DVHMaximumDose,
                                                  'MinDose':dcm.DVHSequence[i].DVHMinimumDose}

                        self.DVHDATA[ROI_name] = dcm.DVHSequence[i].DVHData

                    np.save(self.NEW_DATA_PATH + file_nam + '_DVHSTAT.npy',self.DVHSTAT)
                    np.save(self.NEW_DATA_PATH + file_nam + '_DVHDATA.npy',self.DVHDATA)


                elif dcm.Modality == 'RTPLAN':

                    mu_density = {'Gantry':[],'MLC':[],'JAW':[],'MU':[]}
                    Total_MU = float(dcm.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamMeterset)
                    # MLC = np.zeros([160,len(dcm.BeamSequence[0].ControlPointSequence)])
                    # JAW = np.zeros([2,len(dcm.BeamSequence[0].ControlPointSequence)])
                    print('Total_MU:{}'.format(Total_MU))
                    for j in range(len(dcm.BeamSequence[0].ControlPointSequence)):
                        # extract JAW and MLC
                        mu_density['JAW'].append(np.array(dcm.BeamSequence[0].ControlPointSequence[j].BeamLimitingDevicePositionSequence[0].LeafJawPositions).T)
                        mu_density['MLC'].append(dcm.BeamSequence[0].ControlPointSequence[j].BeamLimitingDevicePositionSequence[1].LeafJawPositions)
                        mu_density['Gantry'].append(dcm.BeamSequence[0].ControlPointSequence[j].GantryAngle)
                        mu_density['MU'].append(Total_MU*float(dcm.BeamSequence[0].ControlPointSequence[j].CumulativeMetersetWeight))
                    # MU = []
                    # for k,ktem in enumerate(MU_):
                    #     if k > 1:
                    #         MU.append((MU_[k]-MU_[k-1])*Total_MU)
                    #     else:
                    #         MU.append(ktem)
                    # print('The length of MU:{}'.format(len(MU)))
                    # with open(self.NEW_DATA_PATH + file_nam + '_Plan.json','w+') as f:
                    #     json.dump(mu_density, f)
                    # mu_density.to_pickle(self.NEW_DATA_PATH + file_nam + '_Plan.pkl')
                    np.save(self.NEW_DATA_PATH + file_nam + '_Plan.npy',mu_density)

            CT = np.array(CT_)
            CT[CT == -2000] = 0

            if slope != 1:
                CT = slope * CT.astype(np.float64)
                CT = img.astype(np.int16)
            CT1 = CT + np.int16(intercept)
            print('The CT array shape: {}'.format(CT1.shape))
            np.save(self.NEW_DATA_PATH + file_nam + '_CT.npy',CT1)

        return mu_density


    def MU_Density_DICOM(self,DICOM_path):
        '''
           Merge all the segment in one fluence map from DICOM via pymedphys tools
           We set the fluence map size is 400 x400 (pixels x pixels), the actual size is 40cm x 40cm.
           So for 1 pixel size is 1 mm.

           The MLC leaf width is 5mm, one bank has 80 leaves.
        '''
        import pymedphys,pydicom
        from pymedphys import Delivery
        # dicom file read
        # DICOM_path = 'C:\\GitFolder\\RL-Application-in-TPS\\AUTO-PLANNING\\AutoTemplateTuning\\projects\\dose prediction\\DATA\\0028\\0028_VMATQATEST.dcm'
        dicom_dataset = pydicom.dcmread(DICOM_path, force=True, stop_before_pixels=True)
        PLAN = Delivery.from_dicom(dicom_dataset)  
        self.PLAN_mu_density = PLAN.mudensity()

        # display MU density
        grid = pymedphys.mudensity.grid()
        pymedphys.mudensity.display(grid,self.PLAN_mu_density)

        return self.PLAN_mu_density

    def MU_Density_TRF(self,TRF_path):
        '''
           Merge all the segment in one fluence map from TRF via pymedphys tool
           We set the fluence map size is 400 x400 (pixels x pixels), the actual size is 40cm x 40cm.
           So for 1 pixel size is 1 mm.

           The MLC leaf width is 5mm, one bank has 80 leaves.
        '''
        import pymedphys
        from pymedphys import Delivery

        # trf file read 
        # trf_file = 'D:\\VM Settings\\003.trf'
        delivery = Delivery.from_logfile(TRF_path)
        self.trf_mu_density = delivery.mudensity() # This is the "fluence"

        # display MU density
        grid = pymedphys.mudensity.grid()
        pymedphys.mudensity.display(grid,self.trf_mu_density)
        
        return self.trf_mu_density

    def display_mu_density_and_difference(self):
        '''
           Display mu difference and plot
        '''
        import pymedphys
        import matplotlib.pyplot as plt
        import numpy as np

        self.grid = pymedphys.mudensity.grid()    # define grid size
        fig = plt.figure(figsize=(5,5))
        # plt.subplot(2,2,1)
        # pymedphys.mudensity.display(grid,self.PLAN_mu_density)
        # plt.title('DICOM')
        # plt.subplot(2,2,2)
        # pymedphys.mudensity.display(grid,self.trf_mu_density)
        # plt.title('TRF files')
        # plt.subplot(2,2,3)
        # pymedphys.mudensity.display_mu_density_diff(grid,self.PLAN_mu_density,self.PLAN_mu_density)

        plt.title('MU difference')
        cmap = "bwr"
        self.diff = self.PLAN_mu_density - self.trf_mu_density
        # if colour_range is None:
        colour_range = np.max(np.abs(self.diff))
        pymedphys.mudensity.display(self.grid,self.diff, 
                                    cmap=cmap,
                                    vmin=-colour_range,
                                    vmax=colour_range)
        print('PLOT finished!')

    def calculate_gamma(self, gamma_options):
        '''
           calculate the gamma of different MU densities
        '''
        import pymedphys
        self.grid = pymedphys.mudensity.grid()    # define grid size

        COORDS = (self.grid["jaw"], self.grid["mlc"])
        gamma = pymedphys.gamma(
        COORDS,
        tuple(map(tuple, self.trf_mu_density)),
        COORDS,
        tuple(map(tuple, self.PLAN_mu_density)),
        **gamma_options)
        
        return gamma


    # def Batch_Generate_MU_density(self,goal_path,
    #                                    TRF_path,
    #                                    DICOM_path):

    #     DICOM_MU_DENSITY = self.MU_Density_DICOM(DICOM_path)
    #     TRF_MU_DENSIT = self.MU_Density_TRF(TRF_path)



    def Beam_Eye_View(self,PLAN):
        '''
           To extract the BEV for each gantry angle.
           
        '''
        return 1

    def DVH_EXTRACTION(self):
        '''
           To extract the DVH data from DICOMs to build up a 
           relationship with anatomy
        '''
    def OVH_calculation(self):
        '''
           To extract strcture set for calculation of overlap volume histogram
        '''







