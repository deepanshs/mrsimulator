# -*- coding: utf-8 -*-
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
from mrsimulator.app.app import app

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


def custom_slider(label="", return_function=None, **kwargs):
    """Custom slider consists of three block -
        1) A label for the slider
        2) Slider bar
        3) A label reflecting the current value of the slider position

        Args:
            label: String with label
            return_function: This function will be applied to the current
                slider value before updating the label.
            kwargs: keyward arguments for dash bootstrap component Input
    """
    id_label = kwargs["id"] + "_label"
    slider = dbc.FormGroup(
        [
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Label(label, color="dark", style={"float": "left"}), width=9
                    ),
                    dbc.Col(dbc.FormText(id=id_label, style={"float": "right"})),
                ]
            ),
            dcc.Slider(**kwargs),
        ]
    )

    @app.callback([Output(id_label, "children")], [Input(kwargs["id"], "value")])
    def update_label(value):
        if return_function is None:
            return [value]
        else:
            return [return_function(value)]

    return slider


def custom_input_group(prepend_label="", append_label="", **kwargs):
    """Custom input group consists of three block -
        1) A prepend label
        2) An Input box
        3) A append label

        Args:
            prepend_label: String to prepend
            append_label: String to append
            kwargs: keyward arguments for dash bootstrap component Input
    """
    return dbc.InputGroup(
        [
            dbc.InputGroupAddon(prepend_label, addon_type="prepend"),
            dbc.Input(pattern="^[-+]?[0-9]*\\.?[0-9]+$", inputMode="numeric", **kwargs),
            dbc.InputGroupAddon(append_label, addon_type="append"),
        ],
        # size="sm",
    )
