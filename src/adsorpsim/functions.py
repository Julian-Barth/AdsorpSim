import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def langmuir_adsorption_carbon(bed_length: float, bed_diameter: float, flow_rate: float, num_segments: int):
    # Physical constants
    R = 8.314  # J/mol·K
    T = 298.15  # K
    P = 101325  # Pa
    initial_molefrac = 400e-6
    initial_conc = initial_molefrac * P / (R * T)  # mol/m³

    # Bed and flow settings
    #bed_length = 1.0  # m
    #bed_diameter = 0.1  # m
    bed_area = np.pi * (bed_diameter / 2) ** 2
    #flow_rate = 0.01  # m³/s
    flow_velocity = flow_rate / bed_area  # m/s

    # Adsorption parameters
    q_max = 2.0  # mol/kg
    K = 1  # 1/(mol/m³)
    k_ads = 0.01  # 1/s
    adsorbent_density = 1000  # kg/m³

    # Discretization
    #num_segments = 100
    dz = bed_length / num_segments

    def initial_conditions():
        C_init = np.zeros(num_segments)
        C_init[0] = initial_conc  # Only inlet has CO₂ at t=0
        q_init = np.zeros(num_segments)
        return np.concatenate([C_init, q_init])


    # Langmuir isotherm
    def q_eq(C):
        return (q_max * K * C) / (1 + K * C)

    # ODE system
    def ode_system(t, y):
        C = y[:num_segments]
        q = y[num_segments:]

        dCdt = np.zeros_like(C)
        dqdt = np.zeros_like(q)

        # Inlet boundary condition
        C_in = initial_conc
        C_upstream = np.concatenate([[C_in], C[:-1]])

        # Advection and adsorption
        dCdz = (C - C_upstream) / dz
        q_eq_vals = q_eq(C)
        dqdt = k_ads * (q_eq_vals - q)
        dCdt = -flow_velocity * dCdz - adsorbent_density * dqdt

        return np.concatenate([dCdt, dqdt])

    # Time span
    t_span = (0, 5000)
    t_eval = np.linspace(*t_span, 1000)

    # Solve
    sol = solve_ivp(
        ode_system,
        t_span,
        initial_conditions(),
        t_eval=t_eval,
        method='BDF',  # good for stiff systems
        vectorized=False,
        rtol=1e-6,
        atol=1e-9
    )

    # Extract outlet CO₂ concentration
    C_sol = sol.y[:num_segments, :]
    outlet_conc = C_sol[-1, :]  # Last segment = outlet

    # Plot breakthrough curve
    plt.figure(figsize=(8, 5))
    plt.plot(sol.t, outlet_conc, label='Outlet CO₂ concentration')
    plt.xlabel('Time (s)')
    plt.ylabel('CO₂ Concentration (mol/m³)')
    plt.title('Breakthrough Curve')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
