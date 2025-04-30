import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import os

from src.adsorpsim.core import Bed, Adsorbent_Langmuir,get_derivative,get_optimal_point  # Adjust if needed

#currentls this is only a "manual mode". The idea would be to have a list of adsorbants we can choose from (at least for Lamgmuir istherms)

# Load adsorbents data
#@st.cache_data()
def download_data():
    "Download the adsorbent database"
    current_file = Path(os.path.abspath(''))
    csv_file = current_file / "data" / "Dataset PPC.csv"
    df = pd.read_csv(csv_file, sep= ";")
    return df

#dowload the data
df = download_data()

# Créer un dictionnaire d'instances, indexé par le nom de l’adsorbant
adsorbants = {}
for _, row in df.iterrows():
    ads = Adsorbent_Langmuir(
        name=row["name"],
        q_max=row["q_max"],
        K0=row["K0"],
        Ea=row["Ea"],
        k_ads=row["k_ads"],
        density=row["density"]
    )
    adsorbants[row["name"]] = ads

#create a list of all the registered adsorbents
list_adsorbents = list(adsorbants.keys())

#Title
st.title('Adsorption Breakthrough Curve Simulator')
st.caption("Practical Proramming In Chemistry - Project")
st.markdown("Julian Barth, Marin Deasgeans, Lucas Eichenbegrer")

# Bed parameters
st.header('Simulation Parameters')
col1, col2 = st.columns([1, 1])
with col1:
    length = st.slider('Bed Length (m)', 0.1, 5.0, 1.0)
    diameter = st.slider('Bed Diameter (m)', 0.01, 0.5, 0.1)
    flow_rate = st.slider('Flow Rate (m³/s)', 0.001, 0.1, 0.01)
with col2:
    num_segments = st.slider('Number of Segments', 10, 500, 100)
    total_time = st.slider('Total Simulation Time (s)', 100, 10000, 3000)
    #définition de l'échelle de la valeure optimale
    optimal_scale = st.slider("Scale of optimal value", 1, 10000, 2000)

# Sidebar inputs
#list of adsorbents  with "Manual Manipulations" at the begining in order to keep the possibility to change parameters manually
st.sidebar.header("List of registered adsorbant")
choix = st.sidebar.selectbox(
    "Choose an adsorbent",
    ["Manual Modifications"]+list_adsorbents,
    index=1  # Par défaut, choisir la première option
)

#Modifs manuelles
if choix=="Manual Modifications":
    # Adsorbent parameters
    st.sidebar.header('Adsorbent Properties')
    q_max = st.sidebar.number_input('q_max (mol/kg)', value=1)
    K0 = st.sidebar.number_input('K (1/(mol/m³))', value=10000,step=10000)
    Ea = st.sidebar.number_input("Ea (J)", value=10000,step=10000)
    k_ads = st.sidebar.number_input('k_ads (1/s)', value=0.01,format="%.4f",step=0.0100)
    density = st.sidebar.number_input('Density (kg/m³)', value=100,step=100)
    # Create adsorbent and bed
    adsorbent = Adsorbent_Langmuir(
    name="Manual adsorbant",
    q_max=q_max,
    K0=K0,
    Ea=Ea,
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
#adsorbents registered in the dataset avec ses propriétés physiques
else:
    #adsorbent parameters
    st.sidebar.header('Adsorbent Properties')
    q_max_fixed = st.sidebar.number_input('q_max (mol/kg)', value=adsorbants[choix].q_max , disabled=True)
    K0_fixed = st.sidebar.number_input('K (1/(mol/m³))', value=adsorbants[choix].K0 , disabled=True,step=10000)
    Ea_fixed = st.sidebar.number_input(("Ea (J)"), value=adsorbants[choix].Ea , disabled=True,step=10000)
    k_ads_fixed = st.sidebar.number_input('k_ads (1/s)', value=adsorbants[choix].k_ads , disabled=True,format="%.4f",step=0.01)
    density_fixed = st.sidebar.number_input('Density (kg/m³)', value=adsorbants[choix].density , disabled=True,step=100)
    # Select adsorbent in the previously created dictionnary and create bed
    adsorbent = adsorbants[choix]
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
optimal_point_x=t[optimal_point]
optimal_point_y=outlet_conc[optimal_point]

# Plotting!!
fig, ax = plt.subplots(figsize=(8, 5))
#Optimal point
ax.plot(optimal_point_x, optimal_point_y, 'rx', markersize=10, label='Optimal point', zorder=2)  # 'r' = rouge, 'x' = croix
#breakthrough curve
ax.plot(t, outlet_conc, label="Outlet CO₂ Concentration",zorder=1)
#template
ax.set_xlabel('Time (s)')
ax.set_ylabel('Outlet CO₂ Concentration (mol/m³)')
ax.set_title('Breakthrough Curve')
ax.legend()
ax.grid(True)
# Toggle to show or not the graph
on_off = st.toggle("Show the graph", value=True)
if on_off:
    st.pyplot(fig)

#test modif