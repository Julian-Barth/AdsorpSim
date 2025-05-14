import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import os
from src.adsorpsim.core import Bed, Adsorbent_Langmuir,download_data,get_percentage_point,add_adsorbent_to_list,plot_the_graph,get_adsorbed_quantity
#the data are loaded: are the data consist of different adsorbents with their physical properties
current_file = Path(os.path.abspath(''))
csv_file = current_file / "data" / "Dataset PPC.csv"
#the data are downloaded
df = download_data(csv_file)
#a dictionnary indexed by the adsorbents names is created, there is an entry for each adsorbent present in the dataset
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
#a list of all the registered adsorbents is created
list_adsorbents = list(adsorbants.keys())



### General inputs:

#Title
st.title('Adsorption Breakthrough Curve Simulator')
st.caption("Practical Proramming In Chemistry - Project")
st.markdown("Julian Barth, Marin Deasgeans, Lucas Eichenbegrer")

#the bed parameters are chosen from the user's interface with defaults settings
st.header('Simulation Parameters')
col1, col2 = st.columns([1, 1])
with col1:
    length = st.slider('Bed Length (m)', 0.1, 5.0, 1.0)
    diameter = st.slider('Bed Diameter (m)', 0.01, 0.5, 0.1)
    flow_rate = st.slider('Flow Rate (m³/s)', 0.001, 0.1, 0.01)
with col2:
    num_segments = st.slider('Number of Segments', 10, 500, 100)
    total_time = st.slider('Total Simulation Time (s)', 100, 10000, 3000)
    T = st.number_input("Temperature (K)", 298.15, step = 10.0)
#a percentage of wanted saturated adsorbent is defined
#this will be further used to calculate the quantity of fixed CO2 depending on the capacity's fraction
percentage_CO2 = st.slider("Percentage of adsorbent saturated in CO₂", 1, 100, 90)



#### Sidebar inputs: 

#the user is invited to choose an adsorbent among the previously prepared list
#the "Manual Manipulations" setting is added to propose the possibility to change parameters manually
st.sidebar.header("List of registered adsorbant")
choix = st.sidebar.selectbox(
    "Choose an adsorbent",
    ["Manual Modifications"]+list_adsorbents,
    index=1
)
#code for the manual modifications setting:
if choix=="Manual Modifications":
    #the adsorbent parameters are chosen from the user's interface with defaults settings
    st.sidebar.header('Adsorbent Properties')
    q_max = st.sidebar.number_input('q_max (mol/kg)', value=1.0, step=1.0)
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
    total_time=total_time,
    adsorbent=adsorbent,
    T=T
)
#code for the adsorbents registered in the dataset
else:
    #adsorbent parameters
    st.sidebar.header('Adsorbent Properties')
    q_max_fixed = st.sidebar.number_input('q_max (mol/kg)', value=adsorbants[choix].q_max , disabled=True, step=1.0)
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
    total_time=total_time,
    adsorbent=adsorbent,
    T=T
)

#Bonus sidebar: the user has the possibility to implement an adsorbents in the dataset directly from the app
st.sidebar.header("Add an adsorbent to the list")
#the constants of the adsorbent has to be listed
with st.sidebar.form("Ajout adsorbent"):
    add_ads_name = st.text_input("name")
    add_ads_q_max = st.number_input("q_max (mol/kg)")
    add_ads_K0 = st.number_input("K0 (1/(mol/m³))")
    add_ads_Ea = st.number_input("Ea (J)")
    add_ads_k_ads = st.number_input("k_ads (1/s)")
    add_ads_density = st.number_input("density (kg/m³)")
    #the use of a form ensures the app will not rerun after each input but only after the button is pressed
    submitted = st.form_submit_button("Add the adsorbent to the list")
#when the button is pressed, the programm will either notify a missing information or add the adsorbent in the dataset
if submitted:
    if add_ads_name=="" or add_ads_q_max==0 or add_ads_K0==0 or add_ads_Ea==0 or add_ads_k_ads==0 or add_ads_density==0:
        st.sidebar.error("there is missing parameter(s)")
    else:
        #add the adsorbent to the .csv file
        add_adsorbent_to_list(csv_file,add_ads_name,add_ads_q_max,add_ads_K0,add_ads_Ea,k_ads,add_ads_density)
        st.sidebar.success("The adsorbent was added to the list, please refresh the page (Ctrl+R)")



#### Plotting:

#different values needed for plotting the graph are calculated using functions that are shown in "core.py"
t, outlet_conc = bed.simulate(plot=False)
pc_point_x, pc_point_y = get_percentage_point(percentage_CO2,t,outlet_conc)
#a toggle to show or not the graph is created
on_off = st.toggle("Show the graph", value=True)
if on_off:
    st.pyplot(plot_the_graph(t,outlet_conc,pc_point_x,pc_point_y))



### Quantities of CO2 adsorbed in one cycle:

st.write("If the experiment is stopped where the red cross is, the following CO₂ quantities would be absorbed")
col11, col22 = st.columns([1, 1])
#in the left part, the quantity of matter of captured CO2 and the acquisition time are shown
with col11:
    #in this tile, the quantity of matter of adsorbed CO2 is calculated and shown
    tile1=st.container(height = 120)
    adsorbed_quantity=get_adsorbed_quantity(t,outlet_conc,pc_point_x,pc_point_y,flow_rate)
    tile1.metric("Quantity of adsorbed CO₂ in one cycle", f"{round(adsorbed_quantity, 2)} [mol]")
    #in this tile, the time of acquisition is shown. It was deducted from the x coordinate of the red cross
    tile3=st.container(height=120)
    tile3.metric("Acquisition time", f"{round(pc_point_x/60)} [min]")
#in the right part, the mass of captured CO2 is shown
with col22:
    tile2=st.container(height = 120)
    tile2.metric("Mass of adsorbed CO₂ in one cycle", f"{round(adsorbed_quantity*0.044009, 2)} [kg]")