import pytest
import numpy as np

from adsorpsim import (
    Adsorbent_Langmuir,
    Bed,
    download_data,
    get_percentage_point,
    plot_the_graph,
    add_adsorbent_to_list,
    get_adsorbed_quantity,
    langmuir_adsorption_carbon
)

@pytest.fixture
def adsorbent():
    return Adsorbent_Langmuir(
        name="TestAdsorbent",
        q_max_CO2=2.0,
        K_CO2=0.5,
        k_ads_CO2=0.01,
        density=1000
    )

@pytest.fixture
def bed(adsorbent):
    return Bed(
        length=1.0,
        diameter=0.1,
        flow_rate=0.01,
        num_segments=10,
        total_time=50,
        adsorbent=adsorbent,
        humidity_percentage=0
    )

def test_adsorbent_repr(adsorbent):
    r = repr(adsorbent)
    assert "TestAdsorbent" in r
    assert "q_max_CO2=2.0" in r

def test_initial_conditions_shape(bed):
    y0 = bed._initial_conditions()
    assert y0.shape[0] == 20  # 10 segments → 10 C_CO2 + 10 q_CO2

def test_simulate_returns_data(bed):
    t, outlet_CO2, outlet_H2O = bed.simulate(plot=False)
    assert isinstance(t, np.ndarray)
    assert isinstance(outlet_CO2, np.ndarray)
    assert outlet_H2O is None
    assert len(t) == bed.total_time

def test_percentage_point_accuracy():
    t = np.linspace(0, 100, 100)
    outlet_conc = np.linspace(0, 0.01624, 100)
    x, y = get_percentage_point(50, t, outlet_conc)
    assert np.isclose(y, 0.00812, atol=1e-4)
    assert 0 <= x <= 100

def test_adsorbed_quantity_computation():
    t = np.linspace(0, 100, 100)
    outlet_conc = 0.01624 - 0.0001 * np.arange(100)  # decreasing
    pc_x, pc_y = 50, 0.01124  # example saturation point
    flow_rate = 0.01
    ads = get_adsorbed_quantity(t, outlet_conc, pc_x, pc_y, flow_rate)
    assert ads > 0