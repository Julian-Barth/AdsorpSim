from adsorpsim.core import Adsorbent_Langmuir, Bed
import numpy as np

def test_adsorbent_eq_loading():
    carbon = Adsorbent_Langmuir(name="TestAC", q_max=2.0, K=1.0, k_ads=0.01, density=1000)

    # Test known values for Langmuir isotherm
    C_test = 1.0  # mol/m³
    q_expected = (2.0 * 1.0 * C_test) / (1 + 1.0 * C_test)  # = 1.0
    q_result = carbon.q_eq(C_test)

    assert np.isclose(q_result, q_expected), f"Expected {q_expected}, got {q_result}"

def test_bed_simulation_runs():
    carbon = Adsorbent_Langmuir(name="TestAC", q_max=2.0, K=1.0, k_ads=0.01, density=1000)
    bed = Bed(length=1.0, diameter=0.1, flow_rate=0.01, num_segments=10, adsorbent=carbon)

    t, C_out = bed.simulate(final_time=100.0)  # Run quick sim

    assert len(t) == len(C_out), "Time and output length mismatch"
    assert np.all(C_out >= 0), "Negative concentration found"
    assert np.all(C_out <= 400), "Concentration exceeds expected range"

def test_initial_conditions_zero_concentration_except_inlet():
    carbon = Adsorbent_Langmuir(name="TestAC", q_max=2.0, K=1.0, k_ads=0.01, density=1000)
    bed = Bed(length=1.0, diameter=0.1, flow_rate=0.01, num_segments=5, adsorbent=carbon)

    C0, _ = bed._initial_conditions()
    assert C0[0] > 0, "Inlet concentration should be non-zero"
    assert np.all(C0[1:] == 0), "Only inlet segment should have non-zero concentration at t=0"