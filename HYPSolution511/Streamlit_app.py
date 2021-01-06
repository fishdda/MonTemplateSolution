import streamlit as st
import pandas as pd
import time

## Set a project title
st.title('Monaco IntelliTemplate V511')


## upload your electronic protocol
uploaded_file = st.file_uploader("Please Upload your Electronic Protocol", type=["xlsx","csv"])

if uploaded_file:
    df = pd.read_excel(uploaded_file)
    st.dataframe(df)

## search if the CT and structure has already been defined




## Machine Learning Tool (PCA+SVR)
uploaded_file2 = st.file_uploader("Please Upload DVH DataBase", type=["xlsx","csv"])
if uploaded_file2:
    DVH_DATA = pd.read_csv("DVH_Parotid_L_Clean.csv")
    st.line_chart(DVH_DATA.iloc[:,1:5])


## Parameters need be entered by users

# dose calculation settings
st.header('Dose Calculation Setting')

dose_alg_option = st.selectbox(
    'Please Select dose calculation algorithm',
    ('Monte Carlo', 'Pencil Beam'))
    
if dose_alg_option == 'Monte Carlo':
    num = st.number_input('Grid Spacing:(cm), Note 0.1-0.8cm')
    st.write('The current number is ',num)

    dose_option = st.selectbox(
        'Statistical Uncertainty(%)',
        ('Per Control Point', 'Per Calculation'))

    if dose_option == 'Per Control Point':
        dose_uncertainty = st.number_input('Per Control Point(%):')
        st.write('The current number is ',dose_uncertainty)
    elif dose_option == 'Per Calculation':
        dose_uncertainty = st.number_input('Per Calculation(%):')
        st.write('The current number is ',dose_uncertainty)


    # sequencing settings
    st.header('Sequencing Setting')

    option = st.selectbox(
        'Select which kind of Delivery Approaches',
        ('Step&Shoot', 'dMLC', 'VMAT'))

    st.write('You selected:', option)

    # VMAT sequencing parameters setting
    if option == "VMAT":

        arc_option = st.selectbox(
            'Maximum Number of Arcs:',
            ('1', '2', '3','4'))

        st.write('You selected:', arc_option)

        number1 = st.number_input('Max.# of Control Points Per Arc:')
        st.write('The current number is ', number1)

        number2 = st.number_input('Min. Segment Width(cm):')
        st.write('The current number is ', number2)  



## generate Monaco Plan Template 
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




## download the Plan Template from Cloud