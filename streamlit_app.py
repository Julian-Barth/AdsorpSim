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
percentage_CO2 = st.slider("Percentage of adsorbent saturated in CO₂", 1, 100, 80)



# Sidebar inputs: adsorbent parameters
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
    q_max = st.sidebar.number_input('q_max (mol/kg)', value=1.0)
    K0 = st.sidebar.number_input('K (1/(mol/m³))', value=10000.0,step=10000.0)
    Ea = st.sidebar.number_input("Ea (J)", value=10000.0,step=10000.0)
    k_ads = st.sidebar.number_input('k_ads (1/s)', value=0.01,format="%.4f",step=0.0100)
    density = st.sidebar.number_input('Density (kg/m³)', value=100.0,step=100.0)
    # Create adsorbent
    adsorbent = Adsorbent_Langmuir(
    name="Manual adsorbant",
    q_max=q_max,
    K0=K0,
    Ea=Ea,
    k_ads=k_ads,
    density=density
    )
    #create bed
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
    K0_fixed = st.sidebar.number_input('K (1/(mol/m³))', value=adsorbants[choix].K0 , disabled=True,step=10000.0)
    Ea_fixed = st.sidebar.number_input(("Ea (J)"), value=adsorbants[choix].Ea , disabled=True,step=10000.0)
    k_ads_fixed = st.sidebar.number_input('k_ads (1/s)', value=adsorbants[choix].k_ads , disabled=True,format="%.4f",step=0.01)
    density_fixed = st.sidebar.number_input('Density (kg/m³)', value=adsorbants[choix].density , disabled=True,step=100.0)
    # Select adsorbent in the previously created dictionnary
    adsorbent = adsorbants[choix]
    #create bed
    bed = Bed(
    length=length,
    diameter=diameter,
    flow_rate=flow_rate,
    num_segments=num_segments,
    adsorbent=adsorbent
)



#values determination /!\ à mettre en fonctions puis mettre les appels dans la partie plotting
# Run simulation
t, outlet_conc = bed.simulate(total_time=total_time, plot=False)
# Calcul de la dérivée
derivative = get_derivative(t, outlet_conc)
#identify the most optimal point
optimal_point = get_optimal_point(derivative,optimal_scale*1e-9)
optimal_point_x=t[optimal_point]
optimal_point_y=outlet_conc[optimal_point]
#identify the point at x percentage of captured carbon
almost_pc_point_y=(percentage_CO2/100*np.max(outlet_conc))
# Trouver l'indice de la valeur la plus proche
idx_proche = np.abs(outlet_conc - almost_pc_point_y).argmin()
# Récupérer la valeur correspondante
pc_point_y = outlet_conc[idx_proche]
index = np.where(outlet_conc == pc_point_y)[0][0]
pc_point_x = t[index]



# Plotting!!
fig, ax = plt.subplots(figsize=(8, 5))
#Optimal point
ax.plot(optimal_point_x, optimal_point_y, 'rx', markersize=10, label='Optimal point', zorder=2)
#percentage point
ax.plot(pc_point_x, pc_point_y, "bx", markersize=10, label="Percentage point",zorder=3)
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



#quantitées adsorbées en un cycle
st.write("If the experiment is stopped where the blue cross is, the following CO₂ quantities would be absorbed")
col11, col22 = st.columns([1, 1])
with col11:
    tile1=st.container(height = 120)
    print(round(pc_point_x))
    sous_outlet_conc=outlet_conc[:round(pc_point_x)+1]
    adsorbed_quantity=0.016291475653314194*pc_point_x*flow_rate-np.sum(sous_outlet_conc)*flow_rate
    tile1.metric("Quantity of adsorbed CO₂ in one cycle", f"{round(adsorbed_quantity, 2)} [mol]")
    tile3=st.container(height=120)
    tile3.metric("Aquisition time", f"{round(pc_point_x/60)} [min]")
with col22:
    tile2=st.container(height = 120)
    tile2.metric("Mass of adsorbed CO₂ in one cycle", f"{round(adsorbed_quantity*44.009, 2)} [g]")



#Bonus sidebar: implementation of adsorbents in the dataset
st.sidebar.header("Add an adsorbent to the list")
with st.sidebar.form("Ajout adsorbent"):
    add_ads_name = st.text_input("name")
    add_ads_q_max = st.number_input("q_max (...)")
    add_ads_K0 = st.number_input("K0 (...)")
    add_ads_Ea = st.number_input("Ea (...)")
    add_ads_k_ads = st.number_input("k_ads (...)")
    add_ads_density = st.number_input("density (...)")

    submitted = st.form_submit_button("Add the adsorbent to the list")

if submitted:
    if add_ads_name=="" or add_ads_q_max==0 or add_ads_K0==0 or add_ads_Ea==0 or add_ads_k_ads==0 or add_ads_density==0:
        st.sidebar.error("a parameter is missing, the adsorbent was not added to the list")
    else:
        st.sidebar.success("The adsorbent was added to the list, please refresh the page (Ctrl+R)")
        #add modification of csv file
        current_file = Path(os.path.abspath(''))
        CSV_PATH = current_file / "data" / "Dataset PPC.csv"
        # Créer un DataFrame avec la nouvelle ligne
        new_row = {
        "name": add_ads_name,
        "q_max": add_ads_q_max,
        "K0": add_ads_K0,
        "Ea": add_ads_Ea,
        "k_ads": add_ads_k_ads,
        "density": add_ads_density
        }

        df = pd.read_csv(CSV_PATH, sep=";")
        df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)

        # Sauvegarder le fichier
        df.to_csv(CSV_PATH, sep=";", index=False)