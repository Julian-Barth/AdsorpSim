import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from pathlib import Path
import pandas as pd
import os

class Adsorbent_Langmuir:
    def __init__(self, name, q_max, K0, Ea, k_ads, density):
        self.name = name
        self.q_max = q_max  # mol/kg
        self.K0 = K0        # pre-exponential factor, 1/(mol/m³)
        self.Ea = Ea        # activation energy in J/mol
        self.k_ads = k_ads  # 1/s
        self.density = density  # kg/m³
        self.R = 8.314      # J/mol·K

    def __repr__(self):
        return f"{self.name} (q_max={self.q_max}, K0={self.K0}, Ea={self.Ea}, k_ads={self.k_ads}, density={self.density})"
    
class Bed: #now we can also initialize a bed with a PREDEFINED adsorbant
    def __init__(self, length, diameter, flow_rate, num_segments, total_time, adsorbent, T):
        self.length = length  # meters
        self.diameter = diameter  # meters
        self.flow_rate = flow_rate  # m³/s
        self.num_segments = num_segments  # discretization segments
        self.total_time = total_time
        self.adsorbent = adsorbent
        self.T = T #Kelvin

        self.area = np.pi * (self.diameter / 2) ** 2
        self.velocity = self.flow_rate / self.area  # m/s
        self.dz = self.length / self.num_segments

        # Inlet conditions, these are specifif to what we do at the CT 
        self.R = 8.314  # J/mol·K
        self.P = 101325  # Pa
        self.initial_molefrac = 400e-6
        self.initial_conc = self.initial_molefrac * self.P / (self.R * self.T)  # mol/m³

    def _initial_conditions(self): #we create a matrix which takes in the concentration of CO2 at any segemnt, in the gas and solid phase. 
        C_init = np.zeros(self.num_segments)
        C_init[0] = self.initial_conc  # Only inlet has CO2 at t=0
        q_init = np.zeros(self.num_segments)
        return np.concatenate([C_init, q_init])

    def _ode_system(self, t, y): 
        # Split the initial conditions back in 2 matrices 
        C = y[:self.num_segments]
        q = y[self.num_segments:]

        dCdt = np.zeros_like(C)
        dqdt = np.zeros_like(q)

        # Inlet boundary condition
        C_in = self.initial_conc
        C_upstream = np.concatenate([[C_in], C[:-1]])

        # Advection and adsorption
        dCdz = (C - C_upstream) / self.dz
        K=self.adsorbent.K0 * np.exp(-self.adsorbent.Ea / (self.adsorbent.R * self.T))
        q_eq_vals=(self.adsorbent.q_max * K * C) / (1 + K * C)
        dqdt = self.adsorbent.k_ads * (q_eq_vals - q)
        dCdt = -self.velocity * dCdz - self.adsorbent.density * dqdt

        return np.concatenate([dCdt, dqdt])

    def simulate(self, plot=False):
        t_span = (0, self.total_time)
        t_eval = np.linspace(*t_span, self.total_time)
        # enbales the integration of the momentum equation 
        sol = solve_ivp(
            self._ode_system, #ode system with as many equations as segments and "t"s
            t_span,
            self._initial_conditions(),
            t_eval=t_eval,
            method='BDF',
            vectorized=False,
            rtol=1e-6,
            atol=1e-9
        )

        C_sol = sol.y[:self.num_segments, :]
        outlet_conc = C_sol[-1, :]

        if plot:
            plt.figure(figsize=(8,5))
            plt.plot(sol.t, outlet_conc, label="Outlet CO₂ Concentration")
            plt.xlabel('Time (s)')
            plt.ylabel('Outlet CO₂ Concentration (mol/m³)')
            plt.title('Breakthrough Curve')
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            plt.show()

        return sol.t, outlet_conc

#the data were not cached as it does not allow the apparition of new asorbent imputted from the app 
def download_data(csv_file):
    "Download the adsorbent database"
    df = pd.read_csv(csv_file, sep= ";")
    return df

def get_percentage_point(percentage:int, t, outlet_conc):
    """
    This function first identifies the outlet CO2 concentration at the desired percentage

    it then finds the index of the nearest value contained in the outlet concentration array 
    (in other words the nearest value taken by the function)

    Once the index is found, we can define and return the point's coordinates
    """
    almost_pc_point_y=(percentage/100*np.max(outlet_conc))

    idx_proche = np.abs(outlet_conc - almost_pc_point_y).argmin()

    pc_point_y = outlet_conc[idx_proche]
    index = np.where(outlet_conc == pc_point_y)[0][0]
    pc_point_x = t[index]
    return pc_point_x, pc_point_y

def plot_the_graph(t,outlet_conc,pc_point_x,pc_point_y):
    """
    This function plots the breakthrough graph
    """
    #the plot is created
    fig, ax = plt.subplots(figsize=(8, 5))
    #the point at wich the desired percentage of adsorbent is saturated in CO2 is plotted
    ax.plot(pc_point_x, pc_point_y, "rx", markersize=10, label="Percentage point",zorder=3)
    #the breakthrough curve is plotted
    ax.plot(t, outlet_conc, label="Outlet CO₂ Concentration",zorder=1)
    #the template of the graph is defined
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Outlet CO₂ Concentration (mol/m³)')
    ax.set_title('Breakthrough Curve')
    ax.legend()
    ax.grid(True)
    return fig

def add_adsorbent_to_list(CSV_PATH, name, q_max, K0, Ea, k_ads, density):
    """
    This function reads a given .csv file 
    and adds a row containing the different needed physical property of an adsorbent
    The physical properties and the name has to be given to the function
    The .csv file is then saved
    """
    #creation of a dataframe with the new line (new adsorbent)
    new_row = {
    "name": name,
    "q_max": q_max,
    "K0": K0,
    "Ea": Ea,
    "k_ads": k_ads,
    "density": density
    }
    #the file is read and the new adsorbent is added
    df = pd.read_csv(CSV_PATH, sep=";")
    df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
    #the file is saved 
    df.to_csv(CSV_PATH, sep=";", index=False)

def get_adsorbed_quantity(t,outlet_conc, pc_point_x, flow_rate):
    """
    
    """
    #the array representing the outlet CO2 concentration is cut to end at the red cross that represent the desired percentage of saturated adsorbent
    #the resulting smaller array's composants are sumed and multiplied by the flowrate to yield the exact quantity of CO2 that went through our system
    #the resulting quantity is substracted from the quantity of CO2 that would have passed through the system if it did not contain any adsorbent
    adsorbed_array=np.array([])
    index = np.where(t == pc_point_x)[0][0]
    for elem in outlet_conc[:round(index)+1]:
        if elem<0.01629:
            adsorbed=0.016301881324379795-elem
            adsorbed_array=np.append(adsorbed_array, adsorbed)
    return sum(adsorbed_array)*flow_rate*pc_point_x

def get_adsorbed_quantity(t,outlet_conc, pc_point_x, flow_rate):
    """
    à revoir
    """
    #the array representing the outlet CO2 concentration is cut to end at the red cross that represent the desired percentage of saturated adsorbent
    #the resulting smaller array's composants are sumed and multiplied by the flowrate to yield the exact quantity of CO2 that went through our system
    #the resulting quantity is substracted from the quantity of CO2 that would have passed through the system if it did not contain any adsorbent
    adsorbed_array=np.array([])
    index = np.where(t == pc_point_x)[0][0]
    for elem in outlet_conc[:round(index)+1]:
        if elem<0.01626:
            adsorbed=0.016301881324379795-elem
            adsorbed_array=np.append(adsorbed_array, adsorbed)
    return sum(adsorbed_array)*flow_rate*pc_point_x




#for testing purposes  
def main():
    carbon = Adsorbent_Langmuir(
        name="Activated Carbon", 
        q_max=2.0, 
        K0=20000,       # example pre-exponential factor
        Ea=25000.0,    # example Ea in J/mol
        k_ads=0.01, 
        density=1000,
       # T=298.15          # custom temperature in Kelvin (optional)
    )
    bed = Bed(
        length=1.0, 
        diameter=0.1, 
        flow_rate=0.01, 
        num_segments=100, 
        total_time=3000,
        adsorbent=carbon,
        T=298.15
    )
    bed.simulate(plot=False)

if __name__ == "__main__":
    main()