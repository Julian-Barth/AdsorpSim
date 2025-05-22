import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import os

from adsorpsim import Bed, Adsorbent_Langmuir
from adsorpsim import download_data,get_percentage_point,add_adsorbent_to_list,plot_the_graph,get_adsorbed_quantity_CO2,get_adsorbed_quantity_H2O,fit_adsorption_parameters_from_df

#the data are loaded: are the data consist of different adsorbents with their physical properties
current_file = Path(os.path.abspath(''))
csv_file = current_file.parents[2] / "AdsorpSim" / "data" / "Adsorbent_data.csv"

#the data are downloaded
df = download_data(csv_file)

#a dictionnary indexed by the adsorbents names is created, there is an entry for each adsorbent present in the dataset
adsorbants = {}
for _, row in df.iterrows():
    ads = Adsorbent_Langmuir(
        name=row["name"],
        q_max_CO2=row["q_max_CO2"],
        K_CO2=row["K_CO2"],
        k_ads_CO2=row["k_ads_CO2"],
        density=row["density"],
        q_max_H2O=row["q_max_H2O"],
        K_H2O=row["K_H2O"],
        k_ads_H2O=row["k_ads_H2O"]
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
    humidity_percentage = st.slider("Percentage of humidity (%)", 0, 100, 0)
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
    index=0
)
st.sidebar.header('Adsorbent Properties')
#code for the manual modifications setting:
if choix=="Manual Modifications":
    #the adsorbent parameters are chosen from the user's interface with defaults settings
    q_max_CO2 = st.sidebar.number_input('Q(max, CO₂) [mol/kg]', value=5.0, step=1.0, key="q_max_CO2")
    K_CO2 = st.sidebar.number_input('K(CO₂) [m³/mol]', value=0.3,format="%.4f", step=0.1, key="K_CO2")
    k_ads_CO2= st.sidebar.number_input('k(ads, H₂O) [1/s]', value=0.01,format="%.4f", step=0.01, key="k_ads_CO2")
    density = st.sidebar.number_input('Density [kg/m³]', value=800.0 ,step=100.0, key="density")
    q_max_H2O = st.sidebar.number_input('Q(max, H₂O) [mol/kg]', value=5.0 , step=1.0, key="q_max_H2O")
    K_H2O= st.sidebar.number_input('K(H₂O) [m³/mol]', value=0.01,format="%.4f", step=0.01, key="K_H2O")
    k_ads_H2O = st.sidebar.number_input('k(ads, H₂O) [1/s]', value=0.1,format="%.4f", step=0.1, key="k_ads_H2O")

    # Create adsorbent
    adsorbent = Adsorbent_Langmuir(
    name="Manual adsorbant",
    q_max_CO2=q_max_CO2,
    K_CO2=K_CO2,
    k_ads_CO2=k_ads_CO2,
    density=density,
    q_max_H2O=q_max_H2O,
    K_H2O=K_H2O,
    k_ads_H2O=k_ads_H2O,
    )
    #create bed
    bed = Bed(
    length=length,
    diameter=diameter,
    flow_rate=flow_rate,
    num_segments=num_segments,
    total_time=total_time,
    adsorbent=adsorbent,
    humidity_percentage=humidity_percentage
)

#code for the adsorbents registered in the dataset
else:
    #adsorbent parameters
    q_max_CO2_ad = st.sidebar.number_input('Q(max, CO₂) [mol/kg]', value=adsorbants[choix].q_max_CO2, step=1.0)
    K_CO2_ad = st.sidebar.number_input('K(CO₂) [m³/mol]', value=adsorbants[choix].K_CO2 , step=0.01,format="%.4f")
    k_ads_CO2_ad = st.sidebar.number_input('k(ads, CO₂) [1/s]', value=adsorbants[choix].k_ads_CO2 ,format="%.4f",step=0.1)
    density_ad = st.sidebar.number_input('Density[kg/m³]', value=adsorbants[choix].density ,step=100.0)
    q_max_H2O_ad = st.sidebar.number_input('Q(max, H₂O) [mol/kg]', value=adsorbants[choix].q_max_H2O , step=1.0)
    K_H2O_ad = st.sidebar.number_input('K (H₂O)[m³/mol]', value=adsorbants[choix].K_H2O, step=0.01,format="%.4f")
    k_ads_H2O_ad = st.sidebar.number_input('k(ads, H₂O) [1/s]', value=adsorbants[choix].k_ads_H2O ,format="%.4f",step=0.1)
    #create modulable adsorbent with the reported parameters
    adsorbent = Adsorbent_Langmuir(
    name="Manual adsorbant",
    q_max_CO2=q_max_CO2_ad,
    K_CO2=K_CO2_ad,
    k_ads_CO2=k_ads_CO2_ad,
    density=density_ad,
    q_max_H2O=q_max_H2O_ad,
    K_H2O=K_H2O_ad,
    k_ads_H2O=k_ads_H2O_ad,
    )
    #create bed
    bed = Bed(
    length=length,
    diameter=diameter,
    flow_rate=flow_rate,
    num_segments=num_segments,
    total_time=total_time,
    adsorbent=adsorbent,
    humidity_percentage=humidity_percentage
)

#Bonus sidebar: the user has the possibility to implement an adsorbents in the dataset directly from the app
st.sidebar.header("Add an adsorbent to the list")
#the constants of the adsorbent has to be listed
with st.sidebar.form("Ajout adsorbent"):
    add_ads_name = st.text_input("name")
    add_ads_q_max_CO2 = st.number_input("Q(max, CO₂) [mol/kg]",step=1.0)
    add_ads_K_CO2 = st.number_input("K(CO₂) [m³/mol]",step=0.1,format="%.4f")
    add_ads_k_ads_CO2 = st.number_input("k(ads, CO₂) [1/s]",step=0.01,format="%.4f")
    add_ads_density = st.number_input("Density [kg/m³]",step=100.0)
    add_ads_q_max_H2O = st.number_input("Q(max, H₂O) [mol/kg] (optional)",step=1.0)
    add_ads_K_H2O = st.number_input("K(H₂O) [m³/mol] (optional)",step=0.1,format="%.4f")
    add_ads_k_ads_H2O = st.number_input("k(ads, H₂O) [1/s] (optional)",step=0.01,format="%.4f")
    #the use of a form ensures the app will not rerun after each input but only after the button is pressed
    submitted = st.form_submit_button("Add the adsorbent to the list")
#when the button is pressed, the programm will either notify a missing information or add the adsorbent in the dataset
if submitted:
    if add_ads_name=="" or add_ads_q_max_CO2==0 or add_ads_K_CO2==0 or add_ads_k_ads_CO2==0 or add_ads_density==0:
        st.sidebar.error("there is missing mandatory parameter(s)")
    elif add_ads_name=="" or add_ads_q_max_CO2<0 or add_ads_K_CO2<0 or add_ads_k_ads_CO2<0 or add_ads_density<0 or add_ads_q_max_H2O<0 or add_ads_K_H2O<0 or add_ads_k_ads_H2O<0:
        st.sidebar.error("the parameters have to be positive")
    elif add_ads_name in list_adsorbents or add_ads_name+" "+"(without H₂O properties)" in list_adsorbents:
        st.sidebar.error("The adsorbent is already in the list")
    elif add_ads_q_max_H2O==0 or add_ads_K_H2O==0 or add_ads_k_ads_H2O==0:
        add_adsorbent_to_list(csv_file,add_ads_name+" "+"(without H₂O properties)",add_ads_q_max_CO2,add_ads_K_CO2,add_ads_k_ads_CO2,add_ads_density)
        st.sidebar.success("The adsorbent was added to the list without specifying its property to adsorb water, please refresh the page (Ctrl+R)")
    else:
    #add the adsorbent to the .csv file
        add_adsorbent_to_list(csv_file,add_ads_name,add_ads_q_max_CO2,add_ads_K_CO2,k_ads_CO2,add_ads_density,add_ads_q_max_H2O,add_ads_K_H2O,k_ads_H2O)
        st.sidebar.success("The adsorbent was added to the list, please refresh the page (Ctrl+R)")



#### Plotting:

#different values needed for plotting the graph are calculated using functions that are shown in "core.py"
t, outlet_CO2, outlet_H2O = bed.simulate()
pc_point_x, pc_point_y = get_percentage_point(percentage_CO2,t,outlet_CO2)
#a toggle to show or not the graph is created
on_off = st.toggle("Show the graph", value=True)
if on_off:
    st.pyplot(plot_the_graph(t,outlet_CO2,outlet_H2O,pc_point_x,pc_point_y))



### Quantities of CO2 adsorbed in one cycle:

st.write("If the experiment is stopped where the red cross is, the following CO₂ quantities would be absorbed")
col11, col22 = st.columns([1, 1])
#in the left part, the quantity of matter of captured CO2 and the acquisition time are shown
with col11:
    #in this tile, the quantity of matter of adsorbed CO2 is calculated and shown
    tile1=st.container(height = 120)
    adsorbed_quantity_CO2=get_adsorbed_quantity_CO2(outlet_CO2,pc_point_x,pc_point_y,flow_rate)
    tile1.metric("Quantity of adsorbed CO₂ in one cycle", f"{round(adsorbed_quantity_CO2, 2)} [mol]")
    #in this tile, the time of acquisition is shown. It was deducted from the x coordinate of the red cross
    tile4=st.container(height=120)
    tile4.metric("Acquisition time", f"{round(pc_point_x/60)} [min]")
with col22:
    tile2=st.container(height = 120)
    tile2.metric("Mass of adsorbed CO₂ in one cycle", f"{round(adsorbed_quantity_CO2*0.044009, 2)} [kg]")
    #in this tile, the mass of adsorbed water is calculated and shown
    tile3=st.container(height = 120)
    tile3.metric("Mass of adsorbed H₂O in one cycle", f"{round(get_adsorbed_quantity_H2O(outlet_CO2,outlet_H2O,humidity_percentage,pc_point_x,pc_point_y,flow_rate)*0.018015, 4)} [kg]")



### adsorbent parameters from csv file:

st.title("Upload your csv file to deduct the adsorbent's parameters")
#show the needed format of the csv file
st.write("The csv file should respect the following format AND use ; as separator")
data = pd.DataFrame(
    [[0, 0.0001],
     [100,0.0003],
     [200,0.0012],
     ["...","..."]],
    columns=["time", "outlet CO₂ concentration"]
)
# Load the csv file
st.dataframe(data, use_container_width=True)
st.markdown("Drag and drop your csv file here, the fitted curve will respect the previously given bed parameters")
st.markdown(" ⚠ Adding a file may significantly increase the running time of the app ⚠")
uploaded_file = st.file_uploader("", type="csv")
#asks the density of the adsorbent to be able to calculate the parameters
if uploaded_file is not None:
    st.write("Please enter the assumed density of the adsorbent.")
    presumed_density=st.number_input("Density [kg/m³]", value=0.0 ,step=100.0, key="presumed density")
    if presumed_density!=0:
        # Calculate the parameters of the adsorbent from the csv file
        df2 = pd.read_csv(uploaded_file,sep=";")
        bed2 = Bed(
            length=length,
            diameter=diameter,
            flow_rate=flow_rate,
            num_segments=num_segments,
            total_time = df2["time"].iloc[-1],
            adsorbent=None,
            humidity_percentage=humidity_percentage
        )
        fitted_adsorbent,fig= fit_adsorption_parameters_from_df(df2,bed2,presumed_density)
        #show the results of the fitted parameters
        st.write("The deducted parameters of the adsorbent are shown below:")
        col1,col2,col3 = st.columns([1, 1, 1])
        with col1:
            tile1=st.container(height = 120)
            tile1.metric("Q(max, CO₂) [mol/kg]", round(fitted_adsorbent.q_max_CO2, 2))
        with col2:
            tile1=st.container(height = 120)
            tile1.metric("K(CO₂) [m³/mol]", round(fitted_adsorbent.K_CO2, 4))
        with col3:
            tile1=st.container(height = 120)
            tile1.metric("k(ads, CO₂) [1/s]", round(fitted_adsorbent.k_ads_CO2, 4))
        #Add the opportunity to directly add the adsorbent to the dataset
        col1,col2 = st.columns([1, 1])
        with col1:
            add_deducted_ads_name = st.text_input("Please name the adsorbent to add it to the list")
        with col2:
            st.write("")
            st.write("")
            if st.button("Add the adsorbent to the list"):
                if add_deducted_ads_name=="":
                    st.sidebar.error("the adsorbent must be named") 
                else:
                    add_adsorbent_to_list(csv_file,add_deducted_ads_name+" "+"(without H₂O properties)",round(fitted_adsorbent.q_max_CO2, 2),round(fitted_adsorbent.K_CO2, 4),round(fitted_adsorbent.k_ads_CO2, 4),presumed_density)
                    st.success("The adsorbent was added to the list without specifying its property to adsorb water, please refresh the page (Ctrl+R)")  
        #show the fitted graph      
        st.pyplot(fig)
        #show the uploaded csv file
        st.write("Aperçu du fichier :")
        st.dataframe(df2)

