import numpy as np
import pandas as pd
import os
from pathlib import Path
from adsorpsim import (
    Adsorbent_Langmuir,
    Bed,
    download_data,
    get_percentage_point,
    get_adsorbed_quantity,
    plot_the_graph,
    add_adsorbent_to_list
)

def test_adsorbent_creation():
    a = Adsorbent_Langmuir("Test", 1.0, 0.1, 0.01, 500)
    assert a.name == "Test"
    assert a.q_max == 1.0
    assert a.K == 0.1
    assert a.k_ads == 0.01
    assert a.density == 500
    assert "Test" in str(a)

def test_bed_initialization():
    a = Adsorbent_Langmuir("Test", 1.0, 0.1, 0.01, 500)
    bed = Bed(1.0, 0.1, 0.01, 10, 100, a, h2o_percentage=0)
    assert bed.length == 1.0
    assert bed.diameter == 0.1
    assert bed.num_segments == 10
    assert bed.initial_conc_CO2 > 0
    assert bed.initial_conc_H2O == 0

def test_simulation_runs_without_h2o():
    a = Adsorbent_Langmuir("Test", 2.0, 1.0, 0.01, 500)
    bed = Bed(1.0, 0.1, 0.01, 20, 1000, a, h2o_percentage=0)
    t, c_out_CO2, c_out_H2O = bed.simulate(plot=False)
    assert len(t) > 0
    assert c_out_CO2[-1] > 0  # Should not be zero
    assert np.allclose(c_out_H2O, 0)

def test_simulation_runs_with_h2o():
    a = Adsorbent_Langmuir("Test", 2.0, 1.0, 0.01, 500)
    bed = Bed(1.0, 0.1, 0.01, 20, 1000, a, h2o_percentage=50)
    t, c_out_CO2, c_out_H2O = bed.simulate(plot=False)
    assert len(t) > 0
    assert np.any(c_out_H2O > 0)

def test_get_percentage_point():
    t = np.linspace(0, 10, 100)
    outlet = np.linspace(0, 0.01624, 100)
    x, y = get_percentage_point(50, t, outlet)
    assert isinstance(x, float)
    assert isinstance(y, float)

def test_get_adsorbed_quantity():
    t = np.linspace(0, 10, 100)
    outlet = np.linspace(0, 0.01624, 100)
    pc_x, pc_y = get_percentage_point(50, t, outlet)
    q = get_adsorbed_quantity(t, outlet, pc_x, pc_y, flow_rate=0.01)
    assert q > 0

def test_add_adsorbent_to_list(tmp_path):
    path = tmp_path / "adsorbents.csv"
    df = pd.DataFrame(columns=["name", "q_max", "K", "k_ads", "density"])
    df.to_csv(path, sep=";", index=False)

    add_adsorbent_to_list(str(path), "TestMaterial", 1.0, 0.5, 0.01, 900)
    df2 = pd.read_csv(path, sep=";")
    assert "TestMaterial" in df2["name"].values

def test_download_data(tmp_path):
    path = tmp_path / "adsorbents.csv"
    df = pd.DataFrame({
        "name": ["Test"],
        "q_max": [1.0],
        "K": [0.1],
        "k_ads": [0.01],
        "density": [500]
    })
    df.to_csv(path, sep=";", index=False)
    loaded_df = download_data(str(path))
    assert not loaded_df.empty
    assert "name" in loaded_df.columns

def test_plot_the_graph():
    t = np.linspace(0, 10, 100)
    conc = np.linspace(0, 0.01624, 100)
    pc_x, pc_y = get_percentage_point(50, t, conc)
    fig = plot_the_graph(t, conc, pc_x, pc_y)
    assert fig is not None