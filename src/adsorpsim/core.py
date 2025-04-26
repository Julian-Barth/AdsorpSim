import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

class Adsorbent_Langmuir:
    def __init__(self, name, q_max, K, k_ads, density):
        self.name = name
        self.q_max = q_max  # mol/kg
        self.K = K          # 1/(mol/m³)
        self.k_ads = k_ads  # 1/s
        self.density = density  # kg/m³

    def equilibrium(self, C):
        return (self.q_max * self.K * C) / (1 + self.K * C)

class Bed:
    def __init__(self, length, diameter, flow_rate, num_segments, adsorbent):
        self.length = length  # meters
        self.diameter = diameter  # meters
        self.flow_rate = flow_rate  # m³/s
        self.num_segments = num_segments  # discretization segments
        self.adsorbent = adsorbent

        self.area = np.pi * (self.diameter / 2) ** 2
        self.velocity = self.flow_rate / self.area  # m/s
        self.dz = self.length / self.num_segments

        # Inlet conditions
        self.R = 8.314  # J/mol·K
        self.T = 298.15  # K
        self.P = 101325  # Pa
        self.initial_molefrac = 400e-6
        self.initial_conc = self.initial_molefrac * self.P / (self.R * self.T)  # mol/m³

    def _initial_conditions(self):
        C_init = np.zeros(self.num_segments)
        C_init[0] = self.initial_conc  # Only inlet has CO₂ at t=0
        q_init = np.zeros(self.num_segments)
        return np.concatenate([C_init, q_init])

    def _ode_system(self, t, y):
        C = y[:self.num_segments]
        q = y[self.num_segments:]

        dCdt = np.zeros_like(C)
        dqdt = np.zeros_like(q)

        # Inlet boundary condition
        C_in = self.initial_conc
        C_upstream = np.concatenate([[C_in], C[:-1]])

        # Advection and adsorption
        dCdz = (C - C_upstream) / self.dz
        q_eq_vals = self.adsorbent.equilibrium(C)
        dqdt = self.adsorbent.k_ads * (q_eq_vals - q)
        dCdt = -self.velocity * dCdz - self.adsorbent.density * dqdt

        return np.concatenate([dCdt, dqdt])

    def simulate(self, total_time=5000, plot=False):
        t_span = (0, total_time)
        t_eval = np.linspace(*t_span, 1000)

        sol = solve_ivp(
            self._ode_system,
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
    
def main():
    carbon = Adsorbent_Langmuir(
        name="Activated Carbon", 
        q_max=2.0, 
        K=1.0, 
        k_ads=0.01, 
        density=1000
    )
    bed = Bed(
        length=1.0, 
        diameter=0.1, 
        flow_rate=0.01, 
        num_segments=100, 
        adsorbent=carbon
    )
    bed.simulate(plot=True)

if __name__ == "__main__":
    main()