# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 18:00:34 2020

@author: xhuae08006
"""
## to run Anatomic_Features' SDAE_FCNN_Model modules ##

from Anatomic_Features import SDAE_FCNN_Model
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd

DVH_DTH_data_path = "E:\\demo\\NBIA_Head_Neck_DATA\\" 
strt_name         = "BrainStem"
epoch             =  500
CNN1D_epoch       =  1000
compress          =  100         # to reduce the original data index spacing of DTH and DVH
file_name_DVH     =  "_DVH_FurtherCleanedUp_ForSDAE.csv"
file_name_DTH     =  "_DTH_FurtherCleanedUp_ForSDAE.csv"
DL_Model = SDAE_FCNN_Model(DVH_DTH_data_path,strt_name)
DL_Model.data_clean_preparation(compress,strt_name+file_name_DVH,strt_name+file_name_DTH)

## DVH SDAE compression 
compressed_data,decoder_data = DL_Model.SDAE_Compression("DVH",epoch)

## Data Preparation for later 1DCNN training
DL_Model._Data_preparation_1DCNN(compress,epoch)

## Predict the test data results 
predict_y,y_test = DL_Model.OneD_CNN_Correlation(strt_name,CNN1D_epoch)

decoder_predict_y = DL_Model.DVH_Reconstruct(strt_name)
        


test_data = DL_Model.DVH_Data[DL_Model.Validat_index_DVH]
decoder_test_y_ = test_data.iloc[::compress,:].T
decoder_test_y_ = pd.DataFrame(decoder_test_y_,index = decoder_predict_y.index,columns = decoder_predict_y.columns)


# summarize history for accuracy
plt.figure()
plt.plot(DL_Model.history.history['accuracy'])
plt.plot(DL_Model.history.history['val_accuracy'])
plt.title('model accuracy')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train', 'validation'], loc='upper left')
plt.grid(b=True)
plt.show()
# summarize history for loss
plt.figure()
plt.plot(DL_Model.history.history['loss'])
plt.plot(DL_Model.history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'validation'], loc='upper left')
plt.grid(b=True)
plt.show()



Test1 = pd.concat([decoder_predict_y.loc[decoder_predict_y.index[0]],
                    decoder_test_y_.loc[decoder_predict_y.index[0]]],axis=1)

Test1.columns = [decoder_predict_y.index[0]+"_SDAE_predicted",decoder_test_y_.index[0]+'_Original_DVH']
Test1['difference'] = np.abs(Test1.iloc[:,0]- Test1.iloc[:,1])

fig = plt.figure()                
Test1.plot()
plt.xlabel('Dose(Gy)')
plt.ylabel('Volume(%)')
plt.title('DVH_'+strt_name)
plt.grid(b=True)
plt.savefig(Test1.columns[0] +pic_name,dpi=200, bbox_inches='tight')
plt.show()




 Test2 = pd.concat([decoder_predict_y.loc[decoder_predict_y.index[1]],
                    decoder_test_y_.loc[decoder_test_y_.index[1]]],axis=1)
 Test2.columns = [decoder_predict_y.index[1]+"_SDAE_predicted",decoder_predict_y.index[1]+'Original_DVH']

 Test2['difference'] = np.abs(Test2.iloc[:,0]- Test2.iloc[:,1])

 fig = plt.figure()                
 Test2.plot()
 plt.xlabel('Dose(Gy)')
 plt.ylabel('Volume(%)')
 plt.title('DVH_ParotidL')
 plt.grid(b=True)
 plt.savefig(Test2.columns[0] +pic_name,dpi=200, bbox_inches='tight')
 plt.show()




 Test3 = pd.concat([decoder_predict_y.loc[decoder_predict_y.index[2]],
                    decoder_test_y_.loc[decoder_test_y.index[2]]],axis=1)
 Test3.columns = [decoder_predict_y.index[2]+"_SDAE_predicted",decoder_predict_y.index[2]+'_Original_DVH']

 Test1['difference'] = np.abs(Test1.iloc[:,0]- Test1.iloc[:,1])
 Test2['difference'] = np.abs(Test2.iloc[:,0]- Test2.iloc[:,1])
 Test3['difference'] = np.abs(Test3.iloc[:,0]- Test3.iloc[:,1])

 pic_name = "_DVH_DecoderVSOriginal_epoch"+str(epoch)+".tif"
  ## check how many information loss after SDAE compression
 fig = plt.figure()                
 Test1.plot()
 plt.xlabel('Dose(Gy)')
 plt.ylabel('Volume(%)')
 plt.title('DVH_ParotidL')
 plt.grid(b=True)
 plt.savefig(Test1.columns[0] +pic_name,dpi=200, bbox_inches='tight')
 plt.show()

 fig = plt.figure()                
 Test2.plot()
 plt.xlabel('Dose(Gy)')
 plt.ylabel('Volume(%)')
 plt.title('DVH_ParotidL')
 plt.grid(b=True)
 plt.savefig(Test2.columns[0] +pic_name,dpi=200, bbox_inches='tight')
 plt.show()

 fig = plt.figure()                
 Test3.plot()
 plt.xlabel('Dose(Gy)')
 plt.ylabel('Volume(%)')
 plt.title('DVH_ParotidL')
 plt.grid(b=True)
 plt.savefig(Test3.columns[0] +pic_name,dpi=200, bbox_inches='tight')
 plt.show()










## Test Data Compressed Results 
# Valid_DVH = DL_Model.Valid_DVH_DATA.iloc[::compress,:].T

# autodecoder_DVH = pd.DataFrame(DL_Model.SC5.inverse_transform(decoder_data),
#                                 index = DL_Model.Valid_DVH_DATA.columns,
#                                 columns = Valid_DVH.columns)

# Test1 = pd.concat([Valid_DVH.loc[Valid_DVH.index[0]],autodecoder_DVH.loc[Valid_DVH.index[0]]],axis=1)
# Test1.columns = [Valid_DVH.index[0]+"_Valid",Valid_DVH.index[0]+'_AutoDecoder']

# Test2 = pd.concat([Valid_DVH.loc[Valid_DVH.index[1]],autodecoder_DVH.loc[Valid_DVH.index[1]]],axis=1)
# Test2.columns = [Valid_DVH.index[1]+"_Valid",Valid_DVH.index[1]+'_AutoDecoder']

# Test3 = pd.concat([Valid_DVH.loc[Valid_DVH.index[2]],autodecoder_DVH.loc[Valid_DVH.index[2]]],axis=1)
# Test3.columns = [Valid_DVH.index[2]+"_Valid",Valid_DVH.index[2]+'_AutoDecoder']

# Test1['difference'] = np.abs(Test1.iloc[:,0]- Test1.iloc[:,1])
# Test2['difference'] = np.abs(Test2.iloc[:,0]- Test2.iloc[:,1])
# Test3['difference'] = np.abs(Test3.iloc[:,0]- Test3.iloc[:,1])

# pic_name = "_DVH_DecoderVSOriginal_epoch"+str(epoch)+".tif"
#  ## check how many information loss after SDAE compression
# fig = plt.figure()                
# Test1.plot()
# plt.xlabel('Dose(Gy)')
# plt.ylabel('Volume(%)')
# plt.title('DVH_ParotidL')
# plt.grid(b=True)
# plt.savefig(Test1.columns[0] +pic_name,dpi=200, bbox_inches='tight')
# plt.show()

# fig = plt.figure()                
# Test2.plot()
# plt.xlabel('Dose(Gy)')
# plt.ylabel('Volume(%)')
# plt.title('DVH_ParotidL')
# plt.grid(b=True)
# plt.savefig(Test2.columns[0] +pic_name,dpi=200, bbox_inches='tight')
# plt.show()

# fig = plt.figure()                
# Test3.plot()
# plt.xlabel('Dose(Gy)')
# plt.ylabel('Volume(%)')
# plt.title('DVH_ParotidL')
# plt.grid(b=True)
# plt.savefig(Test3.columns[0] +pic_name,dpi=200, bbox_inches='tight')
# plt.show()




## DTH SDAE compression 
#compressed_data,decoder_data = DL_Model.SDAE_Compression("DTH",epoch)
#
#Valid_DTH = DL_Model.Valid_DTH_DATA.iloc[::compress,:].T
#
#autodecoder_DTH = pd.DataFrame(DL_Model.SC6.inverse_transform(decoder_data),
#                               index = DL_Model.Valid_DTH_DATA.columns,
#                               columns = Valid_DTH.columns)
#
#Test1 = pd.concat([Valid_DTH.loc[Valid_DTH.index[0]],autodecoder_DTH.loc[Valid_DTH.index[0]]],axis=1)
#Test1.columns = [Valid_DTH.index[0]+"_Valid",Valid_DTH.index[0]+'_AutoDecoder']
#
#Test2 = pd.concat([Valid_DTH.loc[Valid_DTH.index[1]],autodecoder_DTH.loc[Valid_DTH.index[1]]],axis=1)
#Test2.columns = [Valid_DTH.index[1]+"_Valid",Valid_DTH.index[1]+'_AutoDecoder']
#
#Test3 = pd.concat([Valid_DTH.loc[Valid_DTH.index[2]],autodecoder_DTH.loc[Valid_DTH.index[2]]],axis=1)
#Test3.columns = [Valid_DTH.index[2]+"_Valid",Valid_DTH.index[2]+'_AutoDecoder']
#
#Test1[Test1 >=1] = 1
#Test2[Test2 >=1] = 1
#Test3[Test3 >=1] = 1
#
#Test1['difference'] = np.abs(Test1.iloc[:,0]- Test1.iloc[:,1])
#Test2['difference'] = np.abs(Test2.iloc[:,0]- Test2.iloc[:,1])
#Test3['difference'] = np.abs(Test3.iloc[:,0]- Test3.iloc[:,1])
#
#pic_name = "_DTH_DecoderVSOriginal_epoch"+str(epoch)+".tif"
### check how many information loss after SDAE compression
#fig = plt.figure()                
#Test1.plot()
#plt.xlabel('distance(mm)')
#plt.ylabel('Volume(%)')
#plt.title('DTH_ParotidL')
#plt.grid(b=True)
#plt.savefig(Test1.columns[0] +pic_name,dpi=200, bbox_inches='tight')
#plt.show()
#
#fig = plt.figure()                
#Test2.plot()
#plt.xlabel('distance(mm)')
#plt.ylabel('Volume(%)')
#plt.title('DTH_ParotidL')
#plt.grid(b=True)
#plt.savefig(Test2.columns[0] +pic_name,dpi=200, bbox_inches='tight')
#plt.show()
#
#fig = plt.figure()                
#Test3.plot()
#plt.xlabel('distance(mm)')
#plt.ylabel('Volume(%)')
#plt.title('DTH_ParotidL')
#plt.grid(b=True)
#plt.savefig(Test3.columns[0] +pic_name,dpi=200, bbox_inches='tight')
#plt.show()