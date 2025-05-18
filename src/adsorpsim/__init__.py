"""A package to model the adsorption of atmospheric CO₂."""

from __future__ import annotations

from adsorpsim.core import (
    Adsorbent_Langmuir,
    Bed,
    download_data,
    get_percentage_point,
    plot_the_graph,
    add_adsorbent_to_list,
    get_adsorbed_quantity_CO2,
    get_adsorbed_quantity_H2O,
    fit_adsorption_parameters_from_csv,
    load_adsorbent_from_csv
)

__all__ = [
    "Adsorbent_Langmuir",
    "Bed",
    "download_data",
    "get_percentage_point",
    "plot_the_graph",
    "add_adsorbent_to_list",
    "get_adsorbed_quantity_CO2",
    "get_adsorbed_quantity_H2O",
    "fit_adsorption_parameters_from_csv",
    "load_adsorbent_from_csv"
]

__version__ = "0.1.1"
