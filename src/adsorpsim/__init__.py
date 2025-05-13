"""A package to model the adsorption of atmospheric CO2."""

from __future__ import annotations

from .functions import langmuir_adsorption_carbon

from .core import (
    Adsorbent_Langmuir,
    Bed,
    download_data,
    add_adsorbent_to_list,
    get_percentage_point,
    get_adsorbed_quantity,
    plot_the_graph
)

__version__ = "0.0.1"
