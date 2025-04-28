import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

from adsorpsim.core import Bed, Adsorbent_Langmuir  # Adjust if needed

#currentls this is only a "manual mode". The idea would be to have a list of adsorbants we can choose from (at least for Lamgmuir istherms)

st.title('Adsorption Breakthrough Curve Simulator')

# Sidebar inputs
st.sidebar.header('Simulation Parameters')

# Bed parameters
length = st.sidebar.slider('Bed Length (m)', 0.1, 5.0, 1.0)
diameter = st.sidebar.slider('Bed Diameter (m)', 0.01, 0.5, 0.1)
flow_rate = st.sidebar.slider('Flow Rate (m³/s)', 0.001, 0.1, 0.01)
num_segments = st.sidebar.slider('Number of Segments', 10, 500, 100)
total_time = st.sidebar.slider('Total Simulation Time (s)', 100, 10000, 5000)

# Adsorbent parameters
st.sidebar.header('Adsorbent Properties')
q_max = st.sidebar.number_input('q_max (mol/kg)', value=2.0)
K = st.sidebar.number_input('K (1/(mol/m³))', value=1.0)
k_ads = st.sidebar.number_input('k_ads (1/s)', value=0.01)
density = st.sidebar.number_input('Density (kg/m³)', value=1000)

if st.sidebar.button('Run Simulation'):
    # Create adsorbent and bed
    adsorbent = Adsorbent_Langmuir(
        name="Activated Carbon",
        q_max=q_max,
        K=K,
        k_ads=k_ads,
        density=density
    )

    bed = Bed(
        length=length,
        diameter=diameter,
        flow_rate=flow_rate,
        num_segments=num_segments,
        adsorbent=adsorbent
    )

    # Run simulation
    t, outlet_conc = bed.simulate(total_time=total_time, plot=False)

    # Plotting
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t, outlet_conc, label="Outlet CO₂ Concentration")
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Outlet CO₂ Concentration (mol/m³)')
    ax.set_title('Breakthrough Curve')
    ax.legend()
    ax.grid(True)

    st.pyplot(fig)
