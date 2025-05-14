import pytest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Disable GUI backend for plotting
from adsorpsim import (
    Adsorbent_Langmuir, Bed, get_percentage_point,
    plot_the_graph, download_data, add_adsorbent_to_list,
    get_adsorbed_quantity
)

@pytest.fixture
def adsorbent():
    return Adsorbent_Langmuir(
        name="AC_Test",
        q_max_CO2=2.0,
        K_CO2=0.5,
        k_ads_CO2=0.01,
        density=1000,
        q_max_H2O=1.0,
        K_H2O=1,
        k_ads_H2O=0.005
    )

@pytest.fixture
def bed(adsorbent):
    return Bed(
        length=0.5,
        diameter=0.1,
        flow_rate=0.01,
        num_segments=10,
        total_time=100,
        adsorbent=adsorbent,
        humidity_percentage=40
    )

def test_adsorbent_repr(adsorbent):
    text = repr(adsorbent)
    assert "AC_Test" in text
    assert "q_max_CO2=2.0" in text

def test_bed_initial_conditions_with_humidity(bed):
    y0 = bed._initial_conditions()
    assert isinstance(y0, np.ndarray)
    assert len(y0) == 4 * bed.num_segments  # CO2, H2O, q_CO2, q_H2O

def test_bed_simulate(bed):
    t, co2, h2o = bed.simulate()
    assert len(t) == len(co2)
    assert len(t) == len(h2o)
    assert np.all(co2 >= 0)
    assert np.all(h2o >= 0)

def test_get_percentage_point_typical():
    t = np.linspace(0, 100, 100)
    outlet = np.linspace(0, 0.01624, 100)
    x, y = get_percentage_point(50, t, outlet)
    assert isinstance(x, float)
    assert isinstance(y, float)
    assert 0 < x < 100
    assert 0 < y < 0.01624

def test_plot_the_graph_creates_figure():
    t = np.linspace(0, 100, 100)
    outlet = np.linspace(0, 0.01624, 100)
    pc_x, pc_y = get_percentage_point(50, t, outlet)
    fig = plot_the_graph(t, outlet, pc_x, pc_y)
    assert fig is not None
    assert hasattr(fig, "savefig")

def test_get_adsorbed_quantity():
    t = np.linspace(0, 100, 100)
    outlet = np.linspace(0, 0.01624, 100)
    pc_x, pc_y = get_percentage_point(50, t, outlet)
    q_ads = get_adsorbed_quantity(t, outlet, pc_x, pc_y, flow_rate=0.01)
    assert q_ads > 0
    assert isinstance(q_ads, float)

def test_download_data_reads_csv(tmp_path):
    csv_file = tmp_path / "adsorbents.csv"
    csv_file.write_text("name;q_max;K0;Ea;k_ads;density\nAC;2;1;1000;0.01;800\n")
    df = download_data(str(csv_file))
    assert isinstance(df, pd.DataFrame)
    assert "name" in df.columns
    assert df.iloc[0]["name"] == "AC"

def test_add_adsorbent_to_list(tmp_path):
    csv_file = tmp_path / "test_adsorbents.csv"
    df = pd.DataFrame(columns=["name", "q_max", "K0", "Ea", "k_ads", "density"])
    df.to_csv(csv_file, sep=";", index=False)

    add_adsorbent_to_list(
        str(csv_file),
        name="ZeoliteX",
        q_max=1.8,
        K0=0.9,
        Ea=1200,
        k_ads=0.012,
        density=950
    )

    df2 = pd.read_csv(csv_file, sep=";")
    assert "ZeoliteX" in df2["name"].values
    assert df2.shape[0] == 1
