import pytest
import numpy as np
import pandas as pd
import tempfile
from pathlib import Path
from adsorpsim import (
    Adsorbent_Langmuir,
    Bed,
    download_data,
    get_percentage_point,
    plot_the_graph,
    add_adsorbent_to_list,
    get_adsorbed_quantity_CO2,
    get_adsorbed_quantity_H2O,
    fit_adsorption_parameters_from_df,
)

# Fixture: Sample adsorbent for testing
@pytest.fixture
def sample_adsorbent():
    return Adsorbent_Langmuir(
        name="TestAds",
        q_max_CO2=2.0,
        K_CO2=0.5,
        k_ads_CO2=1.0,
        density=1000,
        q_max_H2O=1.0,
        K_H2O=0.1,
        k_ads_H2O=0.5
    )

# Fixture: Sample dry bed using sample adsorbent
@pytest.fixture
def sample_bed(sample_adsorbent):
    return Bed(
        length=1.0,
        diameter=0.1,
        flow_rate=1e-5,
        num_segments=5,
        total_time=10,
        adsorbent=sample_adsorbent
    )

# Test __repr__ method for adsorbent
def test_adsorbent_repr(sample_adsorbent):
    representation = repr(sample_adsorbent)
    assert "TestAds" in representation
    assert "q_max_CO2=2.0" in representation

# Test initial conditions for a dry bed (no humidity)
def test_initial_conditions_dry(sample_bed):
    ic = sample_bed._initial_conditions()
    assert len(ic) == 10  # 5 CO2 concentrations + 5 CO2 adsorbed

# Test initial conditions for humid bed
def test_initial_conditions_humid(sample_adsorbent):
    bed = Bed(1.0, 0.1, 1e-5, 5, 10, sample_adsorbent, humidity_percentage=50)
    ic = bed._initial_conditions()
    assert len(ic) == 20  # CO2 + H2O + their adsorbed quantities

# Test simulation runs without error and returns correct shape
def test_simulate_dry(sample_bed):
    t, outlet_CO2, outlet_H2O = sample_bed.simulate()
    assert len(t) == sample_bed.total_time
    assert len(outlet_CO2) == sample_bed.total_time
    assert outlet_H2O is None

# Test simulation with humidity
def test_simulate_humid(sample_adsorbent):
    bed = Bed(1.0, 0.1, 1e-5, 5, 10, sample_adsorbent, humidity_percentage=50)
    t, outlet_CO2, outlet_H2O = bed.simulate()
    assert outlet_H2O is not None
    assert len(outlet_CO2) == len(outlet_H2O) == bed.total_time

# Test get_percentage_point finds reasonable point
def test_get_percentage_point_valid(sample_bed):
    t, outlet_CO2, _ = sample_bed.simulate()
    x, y = get_percentage_point(50, t, outlet_CO2)
    assert 0 <= x <= max(t)
    assert 0 <= y <= max(outlet_CO2)

# Test graph plotting works
def test_plot_graph_runs(sample_bed):
    t, outlet_CO2, _ = sample_bed.simulate()
    fig = plot_the_graph(t, outlet_CO2, None)
    assert fig is not None

# Test CO2 adsorbed quantity calculation
def test_adsorbed_CO2_quantity(sample_bed):
    t, outlet_CO2, _ = sample_bed.simulate()
    x, y = get_percentage_point(50, t, outlet_CO2)
    q_CO2 = get_adsorbed_quantity_CO2(outlet_CO2, x, y, sample_bed.flow_rate)
    assert q_CO2 >= 0

# Test H2O adsorbed quantity calculation
def test_adsorbed_H2O_quantity(sample_adsorbent):
    bed = Bed(1.0, 0.1, 1e-5, 5, 10, sample_adsorbent, humidity_percentage=50)
    t, outlet_CO2, outlet_H2O = bed.simulate()
    x, y = get_percentage_point(50, t, outlet_CO2)
    q_H2O = get_adsorbed_quantity_H2O(outlet_CO2, outlet_H2O, 50, x, y, bed.flow_rate)
    assert q_H2O >= 0

# Test adding a new adsorbent to a temporary CSV
def test_add_adsorbent_to_csv():
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.csv') as tmpfile:
        df = pd.DataFrame(columns=[
            "name", "q_max_CO2", "K_CO2", "k_ads_CO2", "density",
            "q_max_H2O", "K_H2O", "k_ads_H2O"
        ])
        df.to_csv(tmpfile.name, sep=";", index=False)

        add_adsorbent_to_list(
            tmpfile.name,
            name="NewAds",
            q_max_CO2=1.0,
            K_CO2=0.2,
            k_ads_CO2=0.1,
            density=1000,
            q_max_H2O=0.5,
            K_H2O=0.1,
            k_ads_H2O=0.05
        )

        df_new = pd.read_csv(tmpfile.name, sep=";")
        assert "NewAds" in df_new["name"].values

# Test parameter fitting from dummy experimental data
def test_fit_adsorption_parameters(sample_bed):
    # Create synthetic data based on an existing simulation
    t, outlet_CO2, _ = sample_bed.simulate()
    df_exp = pd.DataFrame({
        "time": t,
        "outlet_CO2": outlet_CO2 + np.random.normal(0, 0.0001, size=len(outlet_CO2))  # add small noise
    })

    fitted_adsorbent, fig = fit_adsorption_parameters_from_df(df_exp, sample_bed)
    assert isinstance(fitted_adsorbent, Adsorbent_Langmuir)
    assert fig is not None
    assert fitted_adsorbent.q_max_CO2 > 0
