import numpy as np
import pytest
from adsorpsim import (
    Adsorbent_Langmuir,
    Bed,
    get_percentage_point,
    get_adsorbed_quantity,
)

@pytest.fixture
def test_adsorbent():
    return Adsorbent_Langmuir(
        name="Test Carbon",
        q_max=2.0,
        K0=20000,
        Ea=25000,
        k_ads=0.01,
        density=1000,
    )

@pytest.fixture
def test_bed(test_adsorbent):
    return Bed(
        length=1.0,
        diameter=0.1,
        flow_rate=0.01,
        num_segments=50,
        total_time=1000,
        adsorbent=test_adsorbent,
        T=298.15
    )

def test_adsorbent_repr(test_adsorbent):
    repr_str = repr(test_adsorbent)
    assert "Test Carbon" in repr_str
    assert "q_max=2.0" in repr_str

def test_bed_simulation_shape(test_bed):
    t, outlet_conc = test_bed.simulate(plot=False)
    assert len(t) == len(outlet_conc)
    assert isinstance(t, np.ndarray)
    assert isinstance(outlet_conc, np.ndarray)
    assert len(outlet_conc) > 0

def test_get_percentage_point_output(test_bed):
    t, outlet_conc = test_bed.simulate(plot=False)
    pc_x, pc_y = get_percentage_point(50, t, outlet_conc)
    assert isinstance(pc_x, (int, float))
    assert isinstance(pc_y, (int, float))
    assert 0 <= pc_y <= max(outlet_conc)

def test_get_adsorbed_quantity(test_bed):
    t, outlet_conc = test_bed.simulate(plot=False)
    pc_x, pc_y = get_percentage_point(50, t, outlet_conc)
    adsorbed = get_adsorbed_quantity(t, outlet_conc, pc_x, pc_y, test_bed.flow_rate)
    assert isinstance(adsorbed, float)
    assert adsorbed >= 0
