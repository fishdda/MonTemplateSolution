import streamlit as st
import pandas as pd
import time
# import pydicom 
import os
from MONV511 import MONV511_UI
from MONV559b import MONV559b_UI
from MONV540 import MONV540_UI
from MONV540_DECODE_V1 import MONV540_UI_DECODE_V1
from MONACO511_MPTD_WEB import MONV511_UI_MPTD

st.sidebar.header('Monaco Plan Template Generation V1.0')

mon_option1 = st.sidebar.radio('Monaco Versions Selection', 
                              ('MONV511','MONV540','MONV559b'))

Mode = st.sidebar.radio('Pareto or Constrained Mode Selection',
                       ('Constrained','Pareto'))

if mon_option1 == 'MONV511':
    MONV511_UI(Mode)
elif mon_option1 == 'MONV559b':
    MONV559b_UI()
elif mon_option1 == 'MONV540':
    MONV540_UI(Mode)


st.sidebar.header('Monaco Plan Template Decode V1.0')

mon_option2 = st.sidebar.selectbox(
    'Monaco Selection',
    ('MONV511','MONV540'))

if mon_option2 == 'MONV540':
    MONV540_UI_DECODE_V1()

st.sidebar.header('Monaco Plan Template Management V1.0')

mon_option3 = st.sidebar.selectbox(
    'Monaco Version',
    ('MONV511','MONV559b'))


st.sidebar.header('Monaco Plan Template Upgrade V1.0 ')

mon_option4 = st.sidebar.selectbox(
    'Monaco TPS Version',
    ('MONV511','MONV559b'))


