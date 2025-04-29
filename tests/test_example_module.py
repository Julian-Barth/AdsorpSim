"""Tests for the AdsorpSim example."""

from adsorpsim import Bed, Adsorbent_Langmuir
import numpy as np

def test_adsorbent_creation():
    carbon = Adsorbent_Langmuir(
        name="TestCarbon",
        q_max=2.0,
        K0=1000.0,
        Ea=10000.0,
        k_ads=0.01,
        density=1000,
        T=298.15
    )
    assert carbon.name == "TestCarbon"
    assert carbon.q_max == 2.0
    assert np.isclose(carbon.K0, 1000.0)
    assert np.isclose(carbon.Ea, 10000.0)
    assert np.isclose(carbon.k_ads, 0.01)
    assert np.isclose(carbon.density, 1000)
    assert isinstance(carbon.K, float)

def test_bed_creation():
    carbon = Adsorbent_Langmuir(
        name="TestCarbon",
        q_max=2.0,
        K0=1000.0,
        Ea=10000.0,
        k_ads=0.01,
        density=1000
    )
    bed = Bed(length=1.0, diameter=0.1, flow_rate=0.01, num_segments=10, adsorbent=carbon)
    assert bed.length == 1.0
    assert bed.diameter == 0.1
    assert bed.flow_rate == 0.01
    assert bed.num_segments == 10
    assert bed.adsorbent.name == "TestCarbon"

def test_bed_simulation_runs():
    carbon = Adsorbent_Langmuir(
        name="TestCarbon",
        q_max=2.0,
        K0=1000.0,
        Ea=10000.0,
        k_ads=0.01,
        density=1000
    )
    bed = Bed(length=1.0, diameter=0.1, flow_rate=0.01, num_segments=10, adsorbent=carbon)
    t, C_out = bed.simulate(plot=False)
    assert len(t) > 0
    assert len(C_out) > 0
    assert t.shape == C_out.shape
