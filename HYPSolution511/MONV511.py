def MONV511_UI(Mode):

    import streamlit as st
    import pandas as pd
    import time
    import pydicom 
    import os
    import base64

    ## Set a project title
    st.title('Monaco IntelliTemplate V511 ('+ Mode +' Mode)')

    ## search if the CT and structure has already been defined



    st.header('Machine Learning Module')
    ## Machine Learning Tool (PCA+SVR)
    st.subheader('PCA+SVR approach')
    uploaded_file2 = st.file_uploader("Please Upload DVH DataBase", type=["csv"])
    if uploaded_file2:
        DVH_DATA = pd.read_csv("DVH_Parotid_L_Clean.csv")
        st.line_chart(DVH_DATA.iloc[:,1:10])

    st.header('E-Protocol Module and Renaming Structure in Monaco')
    ## ====upload your electronic protocol==== ##
    uploaded_file_protocol = st.file_uploader("Please Upload Electronic Protocol", type=["xlsx","csv"])

    if uploaded_file_protocol:
        df = pd.read_excel(uploaded_file_protocol)
        st.dataframe(df)
        patient_id = '00'+str(df.columns[1])     # determine which 

    ## ====upload your rtss.dicom files==== ##
    uploaded_file_rtssDCM = st.file_uploader("Please Upload RT structure DICOM files for renaming(optional)", type=["DCM"])
    if uploaded_file_rtssDCM:
        dcm = pydicom.read_file(uploaded_file_rtssDCM)
        st.text(dcm.Modality)

    
    
    ## ====Parameters need be entered by users==== #

    Parameters_To_Template = {}

    st.header('Monaco Plan Template Generation')

    # Plan Site Checking
    st.subheader('Tumor Sites Checking')
    tumor_option = st.selectbox(
        'Please Checking the Tumor Sites (Currently we only support NPC)',
        ('Prostate','NPC', 'Lung', 'Pancreas', 'Lung SBRT','Liver'))
    st.write('User Select ',tumor_option)

    if tumor_option == 'NPC':

        # Expert or Cancer center
        st.subheader('Expert or Institute Rule')
        rule_option = st.selectbox(
            'Please Select Rule types',
            ('NPC_PekingUnion_6996cGy_Physicist1',
             'NPC_PekingUnion_70Gy_Physicist1', 
             'NPC_PekingUnion_70Gy_Physicist2',
             'NPC_SYSUCC_70Gy_Physicist3',
             'NPC_CAMS_70Gy_Physicist4',
             'NPC_PekingUniversityFirst_70Gy_Physicist4'))
        st.write('User Select ',rule_option)



    # Prescription Checking
    st.subheader('Dose Prescription Setting')

    num_fx = st.number_input('number of fractions:')
    st.write('The current number is ',num_fx)
    Parameters_To_Template[str(num_fx)] = num_fx # send parameters


    fx_dose = st.number_input('Dose per fraction:(Gy)')
    st.write('The current number is ',fx_dose)
    Parameters_To_Template[str(fx_dose)] = fx_dose # send parameters

    # dose calculation settings
    st.subheader('Dose Calculation Setting')

    dose_alg_option = st.selectbox(
        'Please Select dose calculation algorithm',
        ('Monte Carlo', 'Pencil Beam'))
        
    if dose_alg_option == 'Monte Carlo':
        num = st.number_input('Grid Spacing:(cm), Note 0.1-0.8cm')
        st.write('The current number is ',num)
        Parameters_To_Template["dose_grid_spacing"] = num # send parameters
        

        dose_option = st.selectbox('Statistical Uncertainty(%)',
            ('Per Control Point', 'Per Calculation'))

        if dose_option == 'Per Control Point':
            dose_uncertainty = st.number_input('Per Control Point(%):')
            st.write('The current number is ',dose_uncertainty)
            Parameters_To_Template['Per Control Point(%)'] = dose_uncertainty # send parameters
        elif dose_option == 'Per Calculation':
            dose_uncertainty = st.number_input('Per Calculation(%):')
            st.write('The current number is ',dose_uncertainty)
            Parameters_To_Template['Per Calculation(%)'] = dose_uncertainty # send parameters


        # Sequencing settings
        st.subheader('Sequencing Setting')

        option = st.selectbox(
            'Select which kind of Delivery Approaches',
            ('Step&Shoot', 'dMLC', 'VMAT'))

        st.write('You selected:', option)
        Parameters_To_Template['Delivery_Approach'] = option
        delivery_method = Parameters_To_Template['Delivery_Approach']

        Sequencing_Parameters = {}

        # VMAT sequencing parameters setting
        if option == "VMAT":

            Sequencing_Parameters[option] = {}

            # Maximum Number of Arcs in VMAT
            Arc_option = st.selectbox('Maximum Number of Arcs:',
                ('1', '2', '3','4'))
            st.write('You selected:', Arc_option)
            Sequencing_Parameters[option]['Arc Numbers'] = Arc_option # string

            # Maximum Number of Control Points Per Arc
            CP_number = st.number_input('Max.# of Control Points Per Arc:')
            st.write('The current number is ', CP_number)
            Sequencing_Parameters[option]['CP_number per Arc'] = CP_number

            # Minimum Segment Width (cm)
            Min_Seg_Width = st.number_input('Min. Segment Width(cm):')
            st.write('The current number is ', Min_Seg_Width)  
            Sequencing_Parameters[option]['Min. Segment Width(cm)'] = Min_Seg_Width
        
        # dMLC sequencing parameters setting
        elif option == 'dMLC':

            Sequencing_Parameters[option] = {}

            # Maximum Number of Control Points Per Beam
            CP_number_per_beam = st.number_input('Max. # of Control Points Per Beam:')
            st.write('The current number is ', CP_number_per_beam)
            Sequencing_Parameters[option]['CP_number_per_beam'] = CP_number_per_beam # number

            # Minimum Segment Width (cm)
            Min_Seg_Width = st.number_input('Min. Segment Width (cm):')
            st.write('The current number is ', Min_Seg_Width)
            Sequencing_Parameters[option]['Min_Seg_Width'] = Min_Seg_Width

            Max_Sweep_Efficiency = st.checkbox("Max. Sweep Efficiency")    # bool
            Move_Only_Seg        = st.checkbox("Allow Move Only Segments") # bool
            Sequencing_Parameters[option]['Max_Sweep_Efficiency'] = Max_Sweep_Efficiency
            Sequencing_Parameters[option]['Move_Only_Seg'] = Move_Only_Seg

    ## ====Generate Monaco Plan Template==== ## 
    st.header('Click to Batch Generation of Monaco Plan Template')

    if st.button('Start Intelligently Generating a Monaco Template'):
        
        'Starting a long computation...'

        # Add a placeholder
        latest_iteration = st.empty()
        bar = st.progress(0)

        for i in range(100):
            latest_iteration.text(f'TimeCost {i+1}')
            bar.progress(i + 1)
            time.sleep(0.1)

        "...and now we\'re done!"



        ###############################################################################
        #################################   GUI   #####################################


        # delivery_method = Parameters_To_Template['Delivery_Approach']
        # fx = Parameters_To_Template[str(num_fx)]
        # prep_dose=Parameters_To_Template['fx_dose'] * Parameters_To_Template['num_fx']
        # grid_dose = 10* Parameters_To_Template["dose_grid_spacing"]  # mm
        ###############################################################################

        from MONACO511_MPTG_WEB import Initialization_MON511
        # default file path 
        file_path = 'C:/Users/xhuae08006/OneDrive - Elekta/Documents/MonTemplateSolution/HYPSolution511'

        # absolute path for electronic protocol 
        protocol_xlsx_path = os.path.join(file_path,'XH protocol.xlsx')

        # absolute path for structure name changes
        # PT_path = 'C:/Users/Public/Documents/CMS/FocalData/Installation/5~Clinic_XH/1~'
        X = Initialization_MON511(pt_id     = patient_id,
                                  delivery  = delivery_method,
                                  fx        = Parameters_To_Template[str(num_fx)],
                                  prep_dose = Parameters_To_Template[str(fx_dose)] * Parameters_To_Template[str(num_fx)],
                                  grid_dose = 10* Parameters_To_Template["dose_grid_spacing"],
                                  path      = file_path,
                                  protocol_xlsx = protocol_xlsx_path,
                                  sequencing = Sequencing_Parameters)
                
        X.MAIN_GENERATE(tumor_option) # tumor option should be selected by users line #50


    ## download the Plan Template from Cloud
    # st.header('Monaco Plan Templates Download')
    # st.subheader('RTOG Monaco Template Downloads (Cloud Service)')



    st.header("'Monaco Plan Templates Download'")

    text = """\
        There is currently (20191204) no official way of downloading data from Streamlit. See for
        example [Issue 400](https://github.com/streamlit/streamlit/issues/400)

        But I discovered a workaround
        [here](https://github.com/holoviz/panel/issues/839#issuecomment-561538340).

        It's based on the concept of
        [HTML Data URLs](https://developer.mozilla.org/en-US/docs/Web/HTTP/Basics_of_HTTP/Data_URIs)

        You can try it out below for a dataframe csv file download.

        The methodology can be extended to other file types. For inspiration see
        [base64.guru](https://base64.guru/converter/encode/file)
        """
    st.markdown(text)

    data = [(1, 2, 3)]
    # When no file name is given, pandas returns the CSV as a string, nice.
    df = pd.DataFrame(data, columns=["Col1", "Col2", "Col3"])
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
    href = f'<a href="data:file/csv;base64,{b64}">Download CSV File</a> (right-click and save as &lt;some_name&gt;.csv)'
    st.markdown(href, unsafe_allow_html=True)