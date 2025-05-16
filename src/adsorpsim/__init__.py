"""A package to model the adsorption of atmospheric CO2."""

from __future__ import annotations

from .core import (
    Adsorbent_Langmuir,
    Bed,
    download_data,
    get_percentage_point,
    plot_the_graph,
    add_adsorbent_to_list,
    get_adsorbed_quantity,
)

__all__ = [
    "Adsorbent_Langmuir",
    "Bed",
    "download_data",
    "get_percentage_point",
    "plot_the_graph",
    "add_adsorbent_to_list",
    "get_adsorbed_quantity",
    "langmuir_adsorption_carbon",
]

__version__ = "0.0.1"
