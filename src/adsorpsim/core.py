import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from pathlib import Path
import pandas as pd
import os

class Adsorbent_Langmuir:
    def __init__(self, name, q_max_CO2, K_CO2, k_ads_CO2, density, q_max_H2O=0, K_H2O=0, k_ads_H2O=0):
        self.name = name
        self.q_max_CO2 = q_max_CO2
        self.K_CO2 = K_CO2
        self.k_ads_CO2 = k_ads_CO2
        self.density = density
        self.q_max_H2O = q_max_H2O
        self.K_H2O = K_H2O
        self.k_ads_H2O = k_ads_H2O

    def __repr__(self):
        return (f"{self.name} (q_max_CO2={self.q_max_CO2}, K_CO2={self.K_CO2}, k_ads_CO2={self.k_ads_CO2}, "
                f"density={self.density}, q_max_H2O={self.q_max_H2O}, K_H2O={self.K_H2O}, k_ads_H2O={self.k_ads_H2O})")

class Bed:
    def __init__(self, length, diameter, flow_rate, num_segments, total_time, adsorbent, humidity_percentage=0):
        self.length = length
        self.diameter = diameter
        self.flow_rate = flow_rate
        self.num_segments = num_segments
        self.total_time = total_time
        self.adsorbent = adsorbent
        self.humidity_percentage = humidity_percentage  # % humidity (0–100)

        self.area = np.pi * (self.diameter / 2) ** 2
        self.velocity = self.flow_rate / self.area
        self.dz = self.length / self.num_segments

        self.initial_conc_CO2 = 0.01624  # mol/m³

        # Convert relative humidity (%) to water vapor concentration (mol/m³)
        max_H2O_conc = 0.0173  # mol/m³ at 25°C
        self.initial_conc_H2O = (humidity_percentage / 100) * max_H2O_conc

    def _initial_conditions(self):
        C_CO2 = np.zeros(self.num_segments)
        C_CO2[0] = self.initial_conc_CO2
        q_CO2 = np.zeros(self.num_segments)

        if self.initial_conc_H2O > 0:
            C_H2O = np.zeros(self.num_segments)
            C_H2O[0] = self.initial_conc_H2O
            q_H2O = np.zeros(self.num_segments)
            return np.concatenate([C_CO2, C_H2O, q_CO2, q_H2O])
        else:
            return np.concatenate([C_CO2, q_CO2])

    def _ode_system(self, t, y):
        if self.initial_conc_H2O > 0:
            C_CO2 = y[0:self.num_segments]
            C_H2O = y[self.num_segments:2*self.num_segments]
            q_CO2 = y[2*self.num_segments:3*self.num_segments]
            q_H2O = y[3*self.num_segments:]
        else:
            C_CO2 = y[0:self.num_segments]
            q_CO2 = y[self.num_segments:]

        dC_CO2_dt = np.zeros_like(C_CO2)
        C_in_CO2 = self.initial_conc_CO2
        C_up_CO2 = np.concatenate([[C_in_CO2], C_CO2[:-1]])
        dC_CO2_dz = (C_CO2 - C_up_CO2) / self.dz

        K_CO2 = self.adsorbent.K_CO2
        q_eq_CO2 = (self.adsorbent.q_max_CO2 * K_CO2 * C_CO2) / (1 + K_CO2 * C_CO2)
        dq_CO2_dt = self.adsorbent.k_ads_CO2 * (q_eq_CO2 - q_CO2)
        dC_CO2_dt = -self.velocity * dC_CO2_dz - self.adsorbent.density * dq_CO2_dt

        if self.initial_conc_H2O > 0:
            dC_H2O_dt = np.zeros_like(C_H2O)
            C_in_H2O = self.initial_conc_H2O
            C_up_H2O = np.concatenate([[C_in_H2O], C_H2O[:-1]])
            dC_H2O_dz = (C_H2O - C_up_H2O) / self.dz

            K_H2O = self.adsorbent.K_H2O
            q_eq_H2O = (self.adsorbent.q_max_H2O * K_H2O * C_H2O) / (1 + K_H2O * C_H2O)
            dq_H2O_dt = self.adsorbent.k_ads_H2O * (q_eq_H2O - q_H2O)
            dC_H2O_dt = -self.velocity * dC_H2O_dz - self.adsorbent.density * dq_H2O_dt

            return np.concatenate([dC_CO2_dt, dC_H2O_dt, dq_CO2_dt, dq_H2O_dt])
        else:
            return np.concatenate([dC_CO2_dt, dq_CO2_dt])

    def simulate(self, plot=False):
        t_span = (0, self.total_time)
        t_eval = np.linspace(*t_span, self.total_time)

        sol = solve_ivp(
            self._ode_system,
            t_span,
            self._initial_conditions(),
            t_eval=t_eval,
            method='BDF',
            rtol=1e-6,
            atol=1e-9
        )

        C_CO2_sol = sol.y[0:self.num_segments, :]
        outlet_CO2 = C_CO2_sol[-1, :]

        outlet_H2O = None
        if self.initial_conc_H2O > 0:
            C_H2O_sol = sol.y[self.num_segments:2*self.num_segments, :]
            outlet_H2O = C_H2O_sol[-1, :]

        if plot:
            plt.figure(figsize=(8, 5))
            plt.plot(sol.t, outlet_CO2, label="Outlet CO₂")
            if outlet_H2O is not None:
                plt.plot(sol.t, outlet_H2O, label="Outlet H₂O", linestyle='--')
            plt.xlabel('Time (s)')
            plt.ylabel('Outlet Concentration (mol/m³)')
            plt.title('Breakthrough Curve')
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.show()

        return sol.t, outlet_CO2, outlet_H2O

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
    if np.max(outlet_conc)<0.1624:
        almost_pc_point_y=(percentage/100*np.max(outlet_conc))
    else:
        almost_pc_point_y=(percentage/100*0.01624)

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

def get_adsorbed_quantity(t,outlet_conc, pc_point_x, pc_point_y, flow_rate):
    """
    This function will calculate the quantity of adsorbed CO₂ in mol 

    An array containing the concentration of adsorbed CO₂ per m³ is created with all the concentrations from the start to the red cross, (which represents the desired percentage of saturated adsorbent in CO₂)
    
    The component of the array are then summed and multiplied by the acquisition time and the flowrate to retrieve a quantity of matter in moles.
    """
    adsorbed_array=np.array([])
    index = np.where(outlet_conc == pc_point_y)[0][0]
    for elem in outlet_conc[:round(index)+1]:
        if elem<0.01624:
            adsorbed=0.01624-elem
            adsorbed_array=np.append(adsorbed_array, adsorbed)
    return sum(adsorbed_array)*flow_rate*pc_point_x




def main():
    carbon = Adsorbent_Langmuir(
        name="Activated Carbon",
        q_max_CO2=4.5,
        K_CO2=0.2,
        k_ads_CO2=15,
        density=500,
        q_max_H2O=0.15,
        K_H2O=0.05,
        k_ads_H2O=0.005
    )

    bed = Bed(
        length=1,
        diameter=0.15,
        flow_rate=0.01,
        num_segments=100,
        total_time=3000,
        adsorbent=carbon,
        humidity_percentage=50 # Change this to >0 to include humidity
    )

    t, outlet_CO2, outlet_H2O = bed.simulate(plot=True)

if __name__ == "__main__":
    main()