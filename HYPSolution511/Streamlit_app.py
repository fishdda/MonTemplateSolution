import streamlit as st
import pandas as pd
import time
import pydicom 
import os
from MONV511 import MONV511_UI
from MONV559b import MONV559b_UI

st.sidebar.header('MONACO PLAN TEMPLATE GENERATION V1.0')
mon_option = st.sidebar.selectbox(
    'Monaco Versions Selection',
    ('MONV511','MONV551','MONV600','MONV559b'))

Mode = st.sidebar.selectbox(
    'Pareto or Constrained Mode Selection',
    ('Pareto','Constrained'))

if mon_option == 'MONV511':
    MONV511_UI(Mode)
elif mon_option == 'MONV559b':
    MONV559b_UI()


