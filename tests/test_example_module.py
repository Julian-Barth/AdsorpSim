import pytest
import numpy as np
import pandas as pd
from adsorpsim import Adsorbent_Langmuir, Bed, download_data, get_percentage_point, plot_the_graph, add_adsorbent_to_list, get_adsorbed_quantity

# Test Adsorbent_Langmuir class
def test_adsorbent_initialization():
    adsorbent = Adsorbent_Langmuir(
        name="Activated Carbon", 
        q_max=2.0, 
        K0=20000, 
        Ea=25000.0, 
        k_ads=0.01, 
        density=1000
    )
    assert adsorbent.name == "Activated Carbon"
    assert adsorbent.q_max == 2.0
    assert adsorbent.K0 == 20000
    assert adsorbent.Ea == 25000.0
    assert adsorbent.k_ads == 0.01
    assert adsorbent.density == 1000

def test_adsorbent_repr():
    adsorbent = Adsorbent_Langmuir(
        name="Activated Carbon", 
        q_max=2.0, 
        K0=20000, 
        Ea=25000.0, 
        k_ads=0.01, 
        density=1000
    )
    repr_str = repr(adsorbent)
    assert "Activated Carbon" in repr_str
    assert "q_max=2.0" in repr_str
    assert "K0=20000" in repr_str

# Test Bed class
def test_bed_initialization():
    carbon = Adsorbent_Langmuir(name="Activated Carbon", q_max=2.0, K0=20000, Ea=25000.0, k_ads=0.01, density=1000)
    bed = Bed(
        length=1.0, 
        diameter=0.1, 
        flow_rate=0.01, 
        num_segments=100, 
        total_time=3000,
        adsorbent=carbon,
        T=298.15
    )
    assert bed.length == 1.0
    assert bed.diameter == 0.1
    assert bed.flow_rate == 0.01
    assert bed.num_segments == 100
    assert bed.total_time == 3000
    assert bed.adsorbent == carbon
    assert bed.T == 298.15
    assert np.isclose(bed.area, np.pi * (0.1 / 2) ** 2)
    assert np.isclose(bed.velocity, bed.flow_rate / bed.area)
    assert np.isclose(bed.dz, bed.length / bed.num_segments)

def test_bed_initial_conditions():
    carbon = Adsorbent_Langmuir(name="Activated Carbon", q_max=2.0, K0=20000, Ea=25000.0, k_ads=0.01, density=1000)
    bed = Bed(
        length=1.0, 
        diameter=0.1, 
        flow_rate=0.01, 
        num_segments=100, 
        total_time=3000,
        adsorbent=carbon,
        T=298.15
    )
    initial_conditions = bed._initial_conditions()
    assert len(initial_conditions) == 200  # 100 segments for gas and solid
    assert initial_conditions[0] == bed.initial_conc
    assert np.sum(initial_conditions[bed.num_segments:]) == 0

def test_bed_ode_system():
    carbon = Adsorbent_Langmuir(name="Activated Carbon", q_max=2.0, K0=20000, Ea=25000.0, k_ads=0.01, density=1000)
    bed = Bed(
        length=1.0, 
        diameter=0.1, 
        flow_rate=0.01, 
        num_segments=100, 
        total_time=3000,
        adsorbent=carbon,
        T=298.15
    )
    t = 0  # test with a time step
    y = np.zeros(200)  # Initial conditions
    dydt = bed._ode_system(t, y)
    assert dydt.shape == (200,)

# Test Helper Functions
def test_download_data():
    df = download_data("test_adsorbent.csv")
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0

def test_get_percentage_point():
    t = np.linspace(0, 3000, 3000)
    outlet_conc = np.linspace(0, 0.01624, 3000)
    pc_x, pc_y = get_percentage_point(50, t, outlet_conc)
    assert pc_x > 0
    assert pc_y > 0

def test_plot_the_graph():
    t = np.linspace(0, 3000, 3000)
    outlet_conc = np.linspace(0, 0.01624, 3000)
    pc_x, pc_y = get_percentage_point(50, t, outlet_conc)
    fig = plot_the_graph(t, outlet_conc, pc_x, pc_y)
    assert fig is not None

def test_add_adsorbent_to_list():
    # Assuming you have a test CSV file "test_adsorbent.csv"
    add_adsorbent_to_list("test_adsorbent.csv", "New Adsorbent", 2.5, 22000, 28000, 0.02, 1100)
    df = pd.read_csv("test_adsorbent.csv", sep=";")
    assert "New Adsorbent" in df["name"].values

def test_get_adsorbed_quantity():
    t = np.linspace(0, 3000, 3000)
    outlet_conc = np.linspace(0, 0.01624, 3000)
    pc_x, pc_y = get_percentage_point(50, t, outlet_conc)
    adsorbed_quantity = get_adsorbed_quantity(t, outlet_conc, pc_x, pc_y, 0.01)
    assert adsorbed_quantity > 0
