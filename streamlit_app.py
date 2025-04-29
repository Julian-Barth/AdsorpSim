import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

from src.adsorpsim.core import Bed, Adsorbent_Langmuir, get_derivative, get_optimal_point  # Adjust if needed

#currentls this is only a "manual mode". The idea would be to have a list of adsorbants we can choose from (at least for Lamgmuir istherms)

#Title
st.title('How much carbon and how fast can this adsorbant adsorb ??')
st.caption("Practical Proramming In Chemistry - Project")
st.markdown("Julian Barth, Marin Deasgeans, Lucas Eichenbegrer")
st.header('Adsorption Breakthrough Curve Simulator')

# Sidebar inputs
#liste des adsorbants
st.sidebar.header("List of registered adsorbant")
# Choix d'une option unique parmi plusieurs
choix = st.sidebar.selectbox(
    "Select an adsorbant",
    ["Manual Modifications","Activated Carbon", "Option 2", "Option 3"],
    index=0  # Par défaut, choisir la première option
)

# Modifier plusieurs paramètres en fonction des choix
if "Manual Modifications" in choix:
    # Adsorbent parameters
    st.sidebar.header('Adsorbent Properties')
    q_max = st.sidebar.number_input('q_max (mol/kg)', value=2.0)
    K = st.sidebar.number_input('K (1/(mol/m³))', value=1.0)
    k_ads = st.sidebar.number_input('k_ads (1/s)', value=0.01)
    density = st.sidebar.number_input('Density (kg/m³)', value=1000)
    #create the adsorbant
    adsorbent = Adsorbent_Langmuir(
    name="Manual adsorbant",
    q_max=q_max,
    K=K,
    k_ads=k_ads,
    density=density
)
else:
    st.sidebar.header('Adsorbent Properties')
    q_max = st.sidebar.number_input('q_max (mol/kg)', value=2.0 , disabled=True)
    K = st.sidebar.number_input('K (1/(mol/m³))', value=1.0, disabled=True)
    k_ads = st.sidebar.number_input('k_ads (1/s)', value=0.01, disabled=True)
    density = st.sidebar.number_input('Density (kg/m³)', value=1000, disabled=True)
    adsorbent=Adsorbent_Langmuir(
    name = choix
    q_max=q_max,
    K=K,
    k_ads=k_ads,
    density=density
)

# Bed parameters
st.header('Simulation Parameters')
col1, col2 = st.columns([1, 1])
with col1:
    length = st.slider('Bed Length (m)', 0.1, 5.0, 1.0)
    diameter = st.slider('Bed Diameter (m)', 0.01, 0.5, 0.1)
    flow_rate = st.slider('Flow Rate (m³/s)', 0.001, 0.1, 0.01)
with col2:
    num_segments = st.slider('Number of Segments', 10, 500, 100)
    total_time = st.slider('Total Simulation Time (s)', 100, 10000, 5000)
    #définition de l'échelle de la valeure optimale
    optimal_scale = st.slider("Scale of optimal value", 1, 10000, 2000)

# Create adsorbent and bed

bed = Bed(
    length=length,
    diameter=diameter,
    flow_rate=flow_rate,
    num_segments=num_segments,
    adsorbent=adsorbent
)

#values determination
# Run simulation
t, outlet_conc = bed.simulate(total_time=total_time, plot=False)
# Calcul de la dérivée
derivative = get_derivative(t, outlet_conc)
#identify the most optimal point
optimal_point = get_optimal_point(derivative,optimal_scale*1e-9)

# Plotting!!
fig, ax = plt.subplots(figsize=(8, 5))
#Optimal point
ax.plot(t[optimal_point], outlet_conc[optimal_point], 'rx', markersize=10, label='Optimal point', zorder=2)  # 'r' = rouge, 'x' = croix
#breakthrough curve
ax.plot(t, outlet_conc, label="Outlet CO₂ Concentration",zorder=1)
#template
ax.set_xlabel('Time (s)')
ax.set_ylabel('Outlet CO₂ Concentration (mol/m³)')
ax.set_title('Breakthrough Curve')
ax.legend()
ax.grid(True)
st.pyplot(fig)