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
    fit_adsorption_parameters_from_csv,
)
# Create a reusable adsorbent instance for tests
@pytest.fixture
def sample_adsorbent():
    return Adsorbent_Langmuir("TestAds", 2.0, 0.5, 1.0, 1000, 1.0, 0.1, 0.5)

# Create a reusable bed instance using the above adsorbent
@pytest.fixture
def sample_bed(sample_adsorbent):
    return Bed(length=1.0, diameter=0.1, flow_rate=1e-5, num_segments=5, total_time=10, adsorbent=sample_adsorbent)

# Ensure __repr__ for Adsorbent includes the name
def test_adsorbent_repr(sample_adsorbent):
    assert "TestAds" in repr(sample_adsorbent)

# Test initial conditions for a dry bed (no humidity)
def test_bed_initial_conditions_no_humidity(sample_bed):
    ic = sample_bed._initial_conditions()
    assert len(ic) == 10  # Expecting 5 values for CO2 + 5 for adsorbed CO2

# Test initial conditions when humidity is present
def test_bed_initial_conditions_with_humidity(sample_adsorbent):
    bed = Bed(1.0, 0.1, 1e-5, 5, 10, sample_adsorbent, humidity_percentage=50)
    ic = bed._initial_conditions()
    assert len(ic) == 20  # CO2, H2O, and their adsorbed amounts

# Simulate and verify CO2 outlet values for dry conditions
def test_bed_simulate_no_humidity(sample_bed):
    t, outlet_CO2, outlet_H2O = sample_bed.simulate()
    assert len(t) == sample_bed.total_time
    assert outlet_H2O is None  # No H2O expected when humidity is off

# Simulate and verify outlet values with humidity
def test_bed_simulate_with_humidity(sample_adsorbent):
    bed = Bed(1.0, 0.1, 1e-5, 5, 10, sample_adsorbent, humidity_percentage=50)
    t, outlet_CO2, outlet_H2O = bed.simulate()
    assert outlet_H2O is not None  # Expecting H2O outlet data

# Check interpolation works to find a point at a given breakthrough percentage
def test_get_percentage_point():
    t = np.linspace(0, 10, 11)
    outlet = np.linspace(0, 0.01624, 11)
    x, y = get_percentage_point(50, t, outlet)  # Find when outlet reaches 50% of max
    assert isinstance(x, float)
    assert isinstance(y, float)

# Test that plotting function returns a figure (non-null)
def test_plot_the_graph(sample_bed):
    t, outlet_CO2, _ = sample_bed.simulate()
    fig = plot_the_graph(t, outlet_CO2, None, t[3], outlet_CO2[3])
    assert fig is not None

# Test CSV data reading functionality using a temporary file
def test_download_data():
    with tempfile.NamedTemporaryFile(suffix=".csv", mode="w+", delete=False) as f:
        f.write("name;q_max_CO2;K_CO2;k_ads_CO2;density;q_max_H2O;K_H2O;k_ads_H2O\n")
        f.write("Test;1;1;1;1;1;1;1\n")
        f.flush()
        df = download_data(f.name)
    assert not df.empty  # Data should load into a non-empty DataFrame

# Check appending a new adsorbent to file works
def test_add_adsorbent_to_list():
    with tempfile.NamedTemporaryFile(suffix=".csv", mode="w+", delete=False) as f:
        f.write("name;q_max_CO2;K_CO2;k_ads_CO2;density;q_max_H2O;K_H2O;k_ads_H2O\n")
        f.flush()
        add_adsorbent_to_list(f.name, "New", 1, 1, 1, 1, 1, 1, 1)
        df = pd.read_csv(f.name, sep=";")
    assert "New" in df["name"].values

# Calculate CO2 adsorbed from simulation and ensure it’s non-negative
def test_get_adsorbed_quantity_CO2(sample_bed):
    t, outlet_CO2, _ = sample_bed.simulate()
    x, y = get_percentage_point(50, t, outlet_CO2)
    adsorbed = get_adsorbed_quantity_CO2(outlet_CO2, x, y, sample_bed.flow_rate)
    assert adsorbed >= 0

# Calculate H2O adsorbed (when humidity is present) and ensure it’s non-negative
def test_get_adsorbed_quantity_H2O(sample_adsorbent):
    bed = Bed(1.0, 0.1, 1e-5, 5, 10, sample_adsorbent, humidity_percentage=50)
    t, outlet_CO2, outlet_H2O = bed.simulate()
    x, y = get_percentage_point(50, t, outlet_CO2)
    adsorbed = get_adsorbed_quantity_H2O(outlet_CO2, outlet_H2O, 50, x, y, bed.flow_rate)
    assert adsorbed >= 0

# Simulate data and fit parameters to test the fitting function
def test_fit_adsorption_parameters_from_csv(sample_adsorbent):
    bed_template = Bed(1.0, 0.1, 1e-5, 5, 20, sample_adsorbent)
    t_sim, outlet_sim, _ = bed_template.simulate()
    df = pd.DataFrame({'time': t_sim, 'outlet_CO2': outlet_sim})
    fitted_ads, fig = fit_adsorption_parameters_from_csv(df, bed_template)
    assert isinstance(fitted_ads, Adsorbent_Langmuir)  # Should return a fitted model
    assert fig is not None  # Should return a visualization
