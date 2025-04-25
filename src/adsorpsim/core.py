"""Example module to get you started."""


def hello_smiles(smiles: str) -> str:
    """
    Return a greeting string that incorporates the given smiles.

    Parameters
    ----------
    smiles : str
        A text string representing a SMILES (Simplified
        Molecular Input Line Entry System) notation or any string.

    Returns
    -------
    str
        A greeting message incorporating the input smiles.

    Examples
    --------
    >>> hello_smiles("C(=O)O")
    'Hello, C(=O)O!'
    """
    return f"Hello, {smiles}!"


if __name__ == "__main__":
    print(hello_smiles("C(=O)O"))

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


class Adsorbent_Langmuir:
    def __init__(self, name, q_max, K, k_ads, density):
        self.name = name
        self.q_max = q_max          # mol/kg
        self.K = K                  # 1/(mol/m³)
        self.k_ads = k_ads          # 1/s
        self.density = density      # kg/m³

    def q_eq(self, C):
        """Langmuir isotherm"""
        return (self.q_max * self.K * C) / (1 + self.K * C)


class Bed:
    def __init__(self, length, diameter, flow_rate, num_segments, adsorbent):
        self.length = length
        self.diameter = diameter
        self.flow_rate = flow_rate
        self.num_segments = num_segments
        self.adsorbent = adsorbent

        # Constants
        self.R = 8.314     # J/mol·K
        self.T = 298.15    # K
        self.P = 101325    # Pa
        self.initial_molefrac = 400e-6
        self.initial_conc = self.initial_molefrac * self.P / (self.R * self.T)

        # Geometry
        self.bed_area = np.pi * (self.diameter / 2) ** 2
        self.flow_velocity = self.flow_rate / self.bed_area
        self.dz = self.length / self.num_segments

    def initial_conditions(self):
        C_init = np.zeros(self.num_segments)
        C_init[0] = self.initial_conc
        q_init = np.zeros(self.num_segments)
        return np.concatenate([C_init, q_init])

    def ode_system(self, t, y):
        C = y[:self.num_segments]
        q = y[self.num_segments:]

        dCdt = np.zeros_like(C)
        dqdt = np.zeros_like(q)

        # Inlet boundary
        C_in = self.initial_conc
        C_upstream = np.concatenate([[C_in], C[:-1]])

        dCdz = (C - C_upstream) / self.dz
        q_eq_vals = self.adsorbent.q_eq(C)

        dqdt = self.adsorbent.k_ads * (q_eq_vals - q)
        dCdt = -self.flow_velocity * dCdz - self.adsorbent.density * dqdt

        return np.concatenate([dCdt, dqdt])

    def simulate(self, t_final=5000, num_points=1000, plot=True):
        t_span = (0, t_final)
        t_eval = np.linspace(*t_span, num_points)

        sol = solve_ivp(
            self.ode_system,
            t_span,
            self.initial_conditions(),
            t_eval=t_eval,
            method='BDF',
            rtol=1e-6,
            atol=1e-9
        )

        C_sol = sol.y[:self.num_segments, :]
        outlet_conc = C_sol[-1, :]  # Last segment

        if plot:
            plt.figure(figsize=(8, 5))
            plt.plot(sol.t, outlet_conc, label='Outlet CO₂ concentration')
            plt.xlabel('Time (s)')
            plt.ylabel('CO₂ Concentration (mol/m³)')
            plt.title('Breakthrough Curve')
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.show()

        return sol.t, outlet_conc
