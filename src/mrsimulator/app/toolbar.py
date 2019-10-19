# -*- coding: utf-8 -*-
import dash_bootstrap_components as dbc
import dash_html_components as html
from mrsimulator.app.app import app
from dash.dependencies import Input
from dash.dependencies import Output
from dash.dependencies import State
from mrsimulator.app.model_01 import model_01


__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


# tooltips provide additional information.
tooltips = []
tooltip_format = {"placement": "bottom", "delay": {"show": 250, "hide": 10}}
button_format = {"outline": True, "color": "dark", "size": "md"}


def custom_toolbar_buttons(icon="", id=None, tooltip="", **kwargs):
    return dbc.Button(
        [html.I(className=icon), dbc.Tooltip(tooltip, target=id, **tooltip_format)],
        id=id,
        **button_format,
        **kwargs
    )


# Info ------------------------------------------------------------------------------ #
info_button = custom_toolbar_buttons(
    icon="fas fa-info-circle", id="indicator_status", tooltip="Info"
)

# Advance settings ------------------------------------------------------------------ #
setting_button = custom_toolbar_buttons(
    icon="fas fa-cog", id="advance_setting", tooltip="Advance setting"
)

# Show spectrum from individual isotopomers ----------------------------------------- #
decompose_button = custom_toolbar_buttons(
    icon="fas fa-align-left",
    id="decompose",
    tooltip="Show simulation from individual isotopomers",
)


# decompose button callback method
@app.callback(
    [Output("decompose", "active")],
    [Input("decompose", "n_clicks")],
    [State("decompose", "active")],
)
def toggle_decompose_button(n1, status):
    "Toggle decompose button."
    if n1 is not None:
        new_status = True
        if bool(status):
            new_status = False
        return [new_status]
    return [False]


# Button group 1 -------------------------------------------------------------------- #
group_1_buttons = dbc.ButtonGroup(
    [decompose_button, info_button, setting_button, model_01],
    className="btn-group mr-2",
)


# Download dataset ------------------------------------------------------------------ #
# As CSDM
csdm_button = custom_toolbar_buttons(
    icon="fas fa-download", id="download_csdm_button", tooltip="Download file as .csdf"
)


# download buttons callback method
@app.callback(
    [Output("download_csdm_button", "disabled")],
    [Input("nmr_spectrum", "figure")],
    [State("filename_dataset", "children")],
)
def toggle_download_buttons(value, filename_dataset):
    """Toggle download buttons."""
    if filename_dataset in [None, ""]:
        return [True]
    else:
        return [False]


# Button group 1 -------------------------------------------------------------------- #
group2_buttons = dbc.ButtonGroup([html.A(csdm_button, href="", id="download_csdm")])


# toolbar icons --------------------------------------------------------------------- #
toolbar = dbc.Col([group_1_buttons, group2_buttons])
