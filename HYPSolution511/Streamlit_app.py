import streamlit as st
import pandas as pd
import time
import pydicom 
from GT_MONACO511_20200302 import HYP_Editor_MONACO511
import os


## Set a project title
st.title('Monaco IntelliTemplate V511')




## search if the CT and structure has already been defined




## Machine Learning Tool (PCA+SVR)
uploaded_file2 = st.file_uploader("Please Upload DVH DataBase", type=["csv"])
if uploaded_file2:
    DVH_DATA = pd.read_csv("DVH_Parotid_L_Clean.csv")
    st.line_chart(DVH_DATA.iloc[:,1:5])






## ====upload your electronic protocol==== ##
uploaded_file_protocol = st.file_uploader("Please Upload Electronic Protocol", type=["xlsx","csv"])

if uploaded_file_protocol:
    df = pd.read_excel(uploaded_file_protocol)
    st.dataframe(df)

## ====upload your rtss.dicom files==== ##
uploaded_file_rtssDCM = st.file_uploader("Please Upload RT structure DICOM files", type=["DCM"])
if uploaded_file_rtssDCM:
    dcm = pydicom.read_file(uploaded_file_rtssDCM)
    st.text(dcm.Modality)

## ====Parameters need be entered by users==== #

Parameters_To_Template = {}
# Prescription Checking
st.header('Dose Prescription Setting')

num_fx = st.number_input('number of fractions:')
st.write('The current number is ',num_fx)
Parameters_To_Template[str(num_fx)] = num_fx # send parameters


fx_dose = st.number_input('Dose per fraction:(Gy)')
st.write('The current number is ',fx_dose)
Parameters_To_Template[str(fx_dose)] = fx_dose # send parameters

# dose calculation settings
st.header('Dose Calculation Setting')

dose_alg_option = st.selectbox(
    'Please Select dose calculation algorithm',
    ('Monte Carlo', 'Pencil Beam'))
    
if dose_alg_option == 'Monte Carlo':
    num = st.number_input('Grid Spacing:(cm), Note 0.1-0.8cm')
    st.write('The current number is ',num)
    Parameters_To_Template["dose_grid_spacing"] = num # send parameters

    dose_option = st.selectbox(
        'Statistical Uncertainty(%)',
        ('Per Control Point', 'Per Calculation'))

    if dose_option == 'Per Control Point':
        dose_uncertainty = st.number_input('Per Control Point(%):')
        st.write('The current number is ',dose_uncertainty)
        Parameters_To_Template['Per Control Point(%)'] = dose_uncertainty # send parameters
    elif dose_option == 'Per Calculation':
        dose_uncertainty = st.number_input('Per Calculation(%):')
        st.write('The current number is ',dose_uncertainty)
        Parameters_To_Template['Per Calculation(%)'] = dose_uncertainty # send parameters


    # sequencing settings
    st.header('Sequencing Setting')

    option = st.selectbox(
        'Select which kind of Delivery Approaches',
        ('Step&Shoot', 'dMLC', 'VMAT'))

    st.write('You selected:', option)
    Parameters_To_Template['Delivery_Approach'] = option

    # VMAT sequencing parameters setting
    if option == "VMAT":

        arc_option = st.selectbox(
            'Maximum Number of Arcs:',
            ('1', '2', '3','4'))

        st.write('You selected:', arc_option)
        Parameters_To_Template['Arc Numbers'] = arc_option

        CP_number = st.number_input('Max.# of Control Points Per Arc:')
        st.write('The current number is ', CP_number)
        Parameters_To_Template['CP_number per Arc'] = CP_number

        number2 = st.number_input('Min. Segment Width(cm):')
        st.write('The current number is ', number2)  
        Parameters_To_Template['Min. Segment Width(cm)'] = number2



## ====Generate Monaco Plan Template==== ## 
if st.button('Start Intelligently Generating a Monaco Template'):
    
    'Starting a long computation...'

    # Add a placeholder
    latest_iteration = st.empty()
    bar = st.progress(0)

    for i in range(100):
        latest_iteration.text(f'Iteration {i+1}')
        bar.progress(i + 1)
        time.sleep(0.1)

    '...and now we\'re done!(Please Check C:\FocalData\Installation\MonacoTemplate\)'



###############################################################################
#################################   GUI   #####################################
pt_id = '00'+df.columns[1]

delivery_method = Parameters_To_Template['Delivery_Approach']
fx = Parameters_To_Template['num_fx']
prep_dose=Parameters_To_Template['fx_dose'] * Parameters_To_Template['num_fx']
grid_dose = 10* Parameters_To_Template["dose_grid_spacing"]  # mm
###############################################################################


# default file path 
path = 'C:/Users/xhuae08006/OneDrive - Elekta/Documents/MonTemplateSolution/HYPSolution511'

# absolute path for electronic protocol 
protocol_xlsx = os.path.join(path,'XH protocol.xlsx')

# absolute path for structure name changes
PT_path = 'C:/Users/Public/Documents/CMS/FocalData/Installation/5~Clinic_XH/1~'
X = HYP_Editor_MONACO511(pt_id,
                         delivery_method,
                         fx,
                         prep_dose,
                         grid_dose,
                         path,
                         protocol_xlsx,
                         PT_path)
        
X.MAIN_GENERATE('NPC')


## download the Plan Template from Cloud