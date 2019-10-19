# -*- coding: utf-8 -*-
from dash import Dash
import dash_bootstrap_components as dbc
from mrsimulator import Simulator, Dimension
import numpy as np

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


external_scripts = [
    {"src": "https://maxcdn.bootstrapcdn.com/bootstrap/4.1.1/js/bootstrap.min.js"},
    {"src": "https://cdnjs.cloudflare.com/ajax/libs/jquery/3.2.1/jquery.min.js"},
    {
        "src": "https://use.fontawesome.com/releases/v5.0.13/js/solid.js",
        "integrity": (
            "sha384-tzzSw1/Vo+0N5UhStP3bvwWPq+uvzCMfrN1fEFe+xBmv1C/AtVX5K0uZtmcHitFZ"
        ),
        "crossorigin": "anonymous",
    },
    {
        "src": "https://use.fontawesome.com/releases/v5.0.13/js/fontawesome.js",
        "integrity": (
            "sha384-6OIrr52G08NpOFSZdxxz1xdNSndlD4vdcf/q2myIUVO0VsqaGHJsB0RaBE01VTOY"
        ),
        "crossorigin": "anonymous",
    },
    {
        "src": "https://code.jquery.com/jquery-3.3.1.slim.min.js",
        "integrity": (
            "sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo"
        ),
        "crossorigin": "anonymous",
    },
    {
        "src": (
            "https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.0/umd/popper.min.js"
        ),
        "integrity": (
            "sha384-cs/chFZiN24E4KMATLdqdvsezGxaGsi4hLGOzlXwp5UZB1LY//20VyM2taTB4QvJ"
        ),
        "crossorigin": "anonymous",
    },
    {
        "src": "https://stackpath.bootstrapcdn.com/bootstrap/4.1.0/js/bootstrap.min.js",
        "integrity": (
            "sha384-uefMccjFJAIv6A+rW+L4AHf99KvxDjWSu1z9VI8SKNVmz4sk7buKt/6v9KI65qnm"
        ),
        "crossorigin": "anonymous",
    },
]


app = Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    external_scripts=external_scripts,
    meta_tags=[{"name": "viewport", "content": "width=device-width"}],
)
app.config.suppress_callback_exceptions = True

sim = Simulator()
sim.x = np.asarray([-1.0, 1.0])
sim.y = []
sim.original = np.asarray([0.0, 0.0])
sim.nt = 64
sim.spectrum = [Dimension(number_of_points=2, spectral_width=25)]
sim.spectrum_previous = []
sim.new = False
