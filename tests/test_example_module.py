import numpy as np
from adsorpsim.core import Bed, Adsorbent_Langmuir

def test_adsorbent_eq_loading():
    """Test that Langmuir isotherm returns correct equilibrium value."""
    ads = Adsorbent_Langmuir(name="TestAC", q_max=2.0, K=1.0, k_ads=0.01, density=1000)
    C_test = 1.0  # mol/m³
    expected_q_eq = 1.0  # From Langmuir formula: q = (q_max * K * C) / (1 + K * C)
    assert np.isclose(ads.q_eq(C_test), expected_q_eq), "Incorrect Langmuir equilibrium loading"

def test_bed_simulation_runs():
    """Test that the Bed.simulate method runs and returns valid results."""
    carbon = Adsorbent_Langmuir(name="TestAC", q_max=2.0, K=1.0, k_ads=0.01, density=1000)
    bed = Bed(length=1.0, diameter=0.1, flow_rate=0.01, num_segments=10, adsorbent=carbon)
    
    t, C_out = bed.simulate()  # don't pass final_time — use default args

    assert len(t) == len(C_out), "Time and output arrays should be same length"
    assert np.all(C_out >= 0), "All outlet concentrations should be >= 0"
    assert np.max(C_out) < 500, "Unrealistically high concentrations"

def test_simulation_output_monotonic_or_plateau():
    """Basic sanity check: outlet CO₂ shouldn't spike randomly."""
    carbon = Adsorbent_Langmuir(name="TestAC", q_max=2.0, K=1.0, k_ads=0.01, density=1000)
    bed = Bed(length=1.0, diameter=0.1, flow_rate=0.01, num_segments=10, adsorbent=carbon)
    
    t, C_out = bed.simulate()

    # The output should either rise to plateau or stay low (no oscillations)
    diffs = np.diff(C_out)
    num_drops = np.sum(diffs < -1e-3)
    assert num_drops < 3, "Outlet concentration shows unexpected drops"
