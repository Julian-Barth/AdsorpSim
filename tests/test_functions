import pytest
import matplotlib
matplotlib.use("Agg")  # Use non-GUI backend for plotting

import numpy as np
from adsorpsim import langmuir_adsorption_carbon

def test_langmuir_adsorption_carbon_runs_without_error():
    # This test ensures the function completes without exceptions
    langmuir_adsorption_carbon(
        bed_length=0.5,
        bed_diameter=0.1,
        flow_rate=0.01,
        num_segments=10
    )
    assert True  # If it reaches here, it succeeded

def test_simulation_output_values(monkeypatch):
    # Patch plt.show to prevent plot window from popping up
    import matplotlib.pyplot as plt
    monkeypatch.setattr(plt, "show", lambda: None)

    # Run and analyze internal values indirectly
    # Normally, you'd refactor code to expose outputs directly.
    langmuir_adsorption_carbon(
        bed_length=0.2,
        bed_diameter=0.05,
        flow_rate=0.005,
        num_segments=5
    )

    # No exception is a success for now; deep testing requires refactor to expose outputs
    assert True
