import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

class Adsorbent_Langmuir:
    def __init__(self, name, q_max, K0, Ea, k_ads, density, T=298.15):
        self.name = name
        self.q_max = q_max  # mol/kg
        self.K0 = K0        # pre-exponential factor, 1/(mol/m³)
        self.Ea = Ea        # activation energy in J/mol
        self.k_ads = k_ads  # 1/s
        self.density = density  # kg/m³
        self.R = 8.314      # J/mol·K
        self.T = T          # Temperature in K

    @property
    def K(self):
        return self.K0 * np.exp(-self.Ea / (self.R * self.T))

    def equilibrium(self, C):
        K_T = self.K
        return (self.q_max * K_T * C) / (1 + K_T * C)
        
    def __repr__(self):
        return f"{self.name} (q_max={self.q_max}, K={self.K}, k_ads={self.k_ads}, density={self.density})"

class Bed: #now we can also initialize a bed with a PREDEFINED adsorbant
    def __init__(self, length, diameter, flow_rate, num_segments, adsorbent):
        self.length = length  # meters
        self.diameter = diameter  # meters
        self.flow_rate = flow_rate  # m³/s
        self.num_segments = num_segments  # discretization segments
        self.adsorbent = adsorbent

        self.area = np.pi * (self.diameter / 2) ** 2
        self.velocity = self.flow_rate / self.area  # m/s
        self.dz = self.length / self.num_segments

        # Inlet conditions, these are specifif to what we do at the CT 
        self.R = 8.314  # J/mol·K
        self.T = 298.15  # K
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
        q_eq_vals = self.adsorbent.equilibrium(C)
        dqdt = self.adsorbent.k_ads * (q_eq_vals - q)
        dCdt = -self.velocity * dCdz - self.adsorbent.density * dqdt

        return np.concatenate([dCdt, dqdt])

    def simulate(self, total_time=5000, plot=False):
        t_span = (0, total_time)
        t_eval = np.linspace(*t_span, 1000)
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

def get_derivative(t, y):
    """
    Calcule la dérivée numérique dy/dt pour des listes t et y.
    Compatible avec des pas irréguliers.
    """
    t = np.asarray(t)
    y = np.asarray(y)
    dy_dt = np.zeros_like(y)

    # Calcul pour les points internes (centrée)
    dy_dt[1:-1] = (y[2:] - y[:-2]) / (t[2:] - t[:-2])
    
    # Bords : avant et après
    dy_dt[0] = (y[1] - y[0]) / (t[1] - t[0])
    dy_dt[-1] = (y[-1] - y[-2]) / (t[-1] - t[-2])

    return dy_dt

def get_optimal_point(dy_dt, value=2e-6):
    """
    dy_dt : np.array
    value : seuil pour la dérivée

    Retourne l'indice du premier point après le pic de dérivée
    où dy_dt devient inférieur ou égal à value.
    """
    indice_max = np.argmax(dy_dt)
    # Sous-tableau à partir du maximum
    sous_tableau = dy_dt[indice_max:]
    # Trouver les indices où la condition est remplie
    indices = np.where(sous_tableau <= value)[0]
    if len(indices) > 0:
        return indice_max + indices[0]
    else:
        return len(dy_dt)-1  # Aucun point ne correspond
            

#for testing purposes  
def main():
    carbon = Adsorbent_Langmuir(
        name="Activated Carbon", 
        q_max=2.0, 
        K0=20000,       # example pre-exponential factor
        Ea=25000.0,    # example Ea in J/mol
        k_ads=0.01, 
        density=1000,
        T=298.15          # custom temperature in Kelvin (optional)
    )
    bed = Bed(
        length=1.0, 
        diameter=0.1, 
        flow_rate=0.01, 
        num_segments=100, 
        adsorbent=carbon
    )
    bed.simulate(plot=False)

if __name__ == "__main__":
    main()

