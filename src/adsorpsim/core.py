import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from pathlib import Path
import pandas as pd
import os
import warnings

class Adsorbent_Langmuir:
    """
    Represents an adsorbent following Langmuir kinetics.
    """
    def __init__(self, name: str , q_max_CO2 : float, K_CO2: float, k_ads_CO2: float, density: float, q_max_H2O: float=0, K_H2O: float =0, k_ads_H2O: float =0):
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
    """
    Represents a packed bed reactor with discretized segments.
    """
    def __init__(self, length: float, diameter: float, flow_rate: float, num_segments: int, total_time: int, adsorbent : Adsorbent_Langmuir, humidity_percentage: float =0):
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

        if self.initial_conc_H2O != 0:
            C_H2O = np.zeros(self.num_segments)
            C_H2O[0] = self.initial_conc_H2O
            q_H2O = np.zeros(self.num_segments)
            return np.concatenate([C_CO2, C_H2O, q_CO2, q_H2O])
        else:
            return np.concatenate([C_CO2, q_CO2])

    def _ode_system(self, t, y):
        if self.initial_conc_H2O != 0:
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

        if self.initial_conc_H2O != 0:
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

    def simulate(self):
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

        if self.initial_conc_H2O != 0:
            C_H2O_sol = sol.y[self.num_segments:2*self.num_segments, :]
            outlet_H2O = C_H2O_sol[-1, :]
            return sol.t, outlet_CO2, outlet_H2O
        else:
            return sol.t, outlet_CO2, None
        

#the data were not cached as it does not allow the apparition of new asorbent inputted from the app 
def download_data(csv_file):
    "Download the adsorbent database"
    df = pd.read_csv(csv_file, sep= ";")
    return df

def load_adsorbent_from_csv(csv_path, adsorbent_name):
    df = pd.read_csv(csv_path, sep=';')
    row = df[df['name'] == adsorbent_name]

    if row.empty:
        raise ValueError(f"Adsorbent '{adsorbent_name}' not found in {csv_path}.")

    row = row.iloc[0]  # Get the first matching row as Series

    adsorbent = Adsorbent_Langmuir(
        name=row['name'],
        q_max_CO2=row['q_max_CO2'],
        K_CO2=row['K_CO2'],
        k_ads_CO2=row['k_ads_CO2'],
        density=row['density'],
        q_max_H2O=row.get('q_max_H2O', 0),  # default to 0 if not present
        K_H2O=row.get('K_H2O', 0),
        k_ads_H2O=row.get('k_ads_H2O', 0)
    )
    return adsorbent


def get_percentage_point(percentage:float, t, outlet_conc):
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

def plot_the_graph(t,outlet_CO2,outlet_H2O,pc_point_x=None,pc_point_y=None):
    """
    This function plots the breakthrough graph
    """
    #the plot is created
    fig, ax = plt.subplots(figsize=(8, 5))
    if pc_point_x is not None and pc_point_y is not None:
        #the point at wich the desired percentage of adsorbent is saturated in CO2 is plotted
        ax.plot(pc_point_x, pc_point_y, "rx", markersize=10, label="Percentage point",zorder=3)
    #the breakthrough curve is plotted
    ax.plot(t, outlet_CO2, label="Outlet CO₂ Concentration",zorder=1)
    #if the humidity is not null, the outlet concentration of H2O is plotted
    if outlet_H2O is not None:
        ax.plot(t, outlet_H2O, label="Outlet H₂O Concentration", linestyle='--',zorder=2)
    #the template of the graph is defined
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Outlet CO₂ Concentration (mol/m³)')
    ax.set_title('Breakthrough Curve')
    ax.legend()
    ax.grid(True)
    return fig

def add_adsorbent_to_list(CSV_PATH, name, q_max_CO2, K_CO2, k_ads_CO2, density, q_max_H2O=0, K_H2O=0, k_ads_H2O=0):
    """
    This function reads a given .csv file 
    and adds a row containing the different needed physical property of an adsorbent
    The physical properties and the name has to be given to the function
    The .csv file is then saved
    """
    #Verifications for the name
    if name=="":
        raise ValueError("The variable 'name' cannot be empty.")
    if not isinstance(name, str):
        raise TypeError("The variable 'name' must be a string.")
    # Verification for the numeric variables
    numeric_vars = {
        "q_max_CO2": q_max_CO2,
        "K_CO2": K_CO2,
        "k_ads_CO2": k_ads_CO2,
        "density": density,
        "q_max_H2O": q_max_H2O,
        "K_H2O": K_H2O,
        "k_ads_H2O": k_ads_H2O
    }
    for var_name, var_value in numeric_vars.items():
        if not isinstance(var_value, (float, int)):
            raise TypeError(f"The variable '{var_name}' must be a float or an int..")
        if var_value < 0:
            raise ValueError(f"The variable '{var_name}' cannot be negative..")
        if var_name in ["q_max_CO2", "K_CO2", "k_ads_CO2", "density"] and var_value == 0:
            raise ValueError(f"The variable '{var_name}' cannot be zero.")
    #the file is read 
    df = pd.read_csv(CSV_PATH, sep=";")
    # Check if the adsorbent already exists
    if name in df['name'].values:
        raise ValueError(f"The adsorbent '{name}' is already in the database.")
    #creation of a dataframe with the new line (new adsorbent)
    new_row = {
    "name": name,
    "q_max_CO2": q_max_CO2,
    "K_CO2": K_CO2,
    "k_ads_CO2": k_ads_CO2,
    "density": density,
    "q_max_H2O": q_max_H2O,
    "K_H2O": K_H2O,
    "k_ads_H2O": k_ads_H2O
    }
    # the new adsorbent is added
    df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
    #the file is saved 
    df.to_csv(CSV_PATH, sep=";", index=False)

def get_adsorbed_quantity_CO2(outlet_conc, pc_point_x, pc_point_y, flow_rate):
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

def get_adsorbed_quantity_H2O(outlet_CO2,outlet_H2O, humidity_precentage, pc_point_x, pc_point_y, flow_rate):
    """
    This function will calculate the quantity of adsorbed H₂O in mol 

    An array containing the concentration of adsorbed H₂O per m³ is created with all the concentrations from the start to the red cross, (which represents the desired percentage of saturated adsorbent in CO₂)
    
    The component of the array are then summed and multiplied by the acquisition time and the flowrate to retrieve a quantity of matter in moles.
    """
    if outlet_H2O is not None:
        adsorbed_array=np.array([])
        index = np.where(outlet_CO2 == pc_point_y)[0][0]
        max_value = 0.0173 * (humidity_precentage / 100)  # max H2O concentration at given humidity
        for elem in outlet_H2O[:round(index)+1]:
            if elem<max_value:
                adsorbed=max_value-elem
                adsorbed_array=np.append(adsorbed_array, adsorbed)
        return sum(adsorbed_array)*flow_rate*pc_point_x
    else:
        return 0

def fit_adsorption_parameters_from_df(df, bed_template,  assumed_density = None, initial_guess=[4.0, 0.2, 1]):
    """
    Fit the Langmuir adsorption parameters to experimental CO2 breakthrough data.

    Parameters:
        df : dataframe containing 'time' and 'outlet_CO2' columns.
        bed_template (Bed): A Bed object with all fixed parameters except the adsorbent (can use dummy adsorbent initially).
        initial_guess (list): [q_max_CO2, K_CO2, k_ads_CO2]

    Returns:
        fitted_adsorbent (Adsorbent_Langmuir): Fitted adsorbent object.
        fig (matplotlib figure): Figure of the fitted breakthrough curve.
    """
    
    if bed_template.adsorbent == None:

        bed_template.adsorbent = Adsorbent_Langmuir(
        name="Dummy Adsorbant",
        q_max_CO2=5.0,
        K_CO2=0.2,
        k_ads_CO2=0.02,
        density=assumed_density
)
    # Load experimental data
    t_exp = df['time'].values
    outlet_CO2_exp = df['outlet_CO2'].values

    def loss(params):
        q_max, K, k_ads = params

        ads = Adsorbent_Langmuir(
            name="Fitted",
            q_max_CO2=q_max,
            K_CO2=K,
            k_ads_CO2=k_ads,
            density=bed_template.adsorbent.density,
        )
        bed = Bed(
            length=bed_template.length,
            diameter=bed_template.diameter,
            flow_rate=bed_template.flow_rate,
            num_segments=bed_template.num_segments,
            total_time=bed_template.total_time,
            adsorbent=ads,
            humidity_percentage=bed_template.humidity_percentage
        )

        try:
            t_model, outlet_model, _ = bed.simulate()
            outlet_interp = np.interp(t_exp, t_model, outlet_model)
            error = np.mean((outlet_interp - outlet_CO2_exp)**2)
            return error
        except Exception as e:
            return 1e6  # penalize failed simulations

    # Optimization
    result = minimize(loss, initial_guess, method='Nelder-Mead')

    q_max_opt, K_opt, k_ads_opt = result.x

    fitted_adsorbent = Adsorbent_Langmuir(
        name="Fitted_Adsorbent",
        q_max_CO2=q_max_opt,
        K_CO2=K_opt,
        k_ads_CO2=k_ads_opt,
        density=bed_template.adsorbent.density
    )

    # Plot experimental vs model
    bed_template.adsorbent = fitted_adsorbent
    t_sim, outlet_sim, _ = bed_template.simulate()
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t_sim, outlet_sim, label="Fitted model",zorder=1)
    ax.scatter(t_exp, outlet_CO2_exp, label="Experimental points",color="red", marker="o",zorder=2)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Outlet CO₂ Concentration (mol/m³)')
    ax.set_title('Breakthrough Curve')
    ax.legend()
    ax.grid(True)
    
    return fitted_adsorbent,fig