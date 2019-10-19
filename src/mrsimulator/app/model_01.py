# -*- coding: utf-8 -*-
"""Model page layout and callbacks"""

import dash_bootstrap_components as dbc
import dash_core_components as dcc
from mrsimulator.app.app import app
from dash.dependencies import Input
from dash.dependencies import Output

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


# Line 1 ----------------------------------------------------------------------
# number of orientation used in averaging
model_line_1 = dbc.Row(
    [
        dbc.Col(dbc.Label("Integration density")),
        dbc.Col(
            dbc.Input(
                type="number",
                value=70,
                # type="slider",
                min=0,
                max=4096,
                step=1,
                id="averaging_quality",
            )
        ),
    ]
)
model_line_1_info = dbc.Label(
    size="sm", id="number_of_averaging_points", style={"color": "#566573"}
)


@app.callback(
    Output("number_of_averaging_points", "children"),
    [Input("averaging_quality", "value")],
)
def update_number_of_orientations(value):
    """
    Update the number of orientation for powder averaging.
    Option for advance modal.
    """
    ori = 2 * (value + 1) * (value + 2)
    return f"Averaging over {ori} orientations.".format(value)


# Line 2 ----------------------------------------------------------------------
# integration volume. Options are Octant, Hemisphere, Sphere
model_line_2 = dbc.Row(
    [
        dbc.Col(dbc.Label("Integration volume")),
        dbc.Col(
            dcc.Dropdown(
                id="n_octants",
                options=[
                    {"label": "Octant", "value": 0},
                    {"label": "Hemisphere", "value": 1},
                    {"label": "Sphere", "value": 2},
                ],
                value=1,
                clearable=False,
            )
        ),
    ]
)


# Layout ----------------------------------------------------------------------
# model user-interface
model_01 = dbc.Modal(
    [
        dbc.ModalHeader("Advance setting"),
        dbc.ModalBody(dbc.FormGroup([model_line_1, model_line_1_info, model_line_2])),
        dbc.ModalFooter(
            dbc.Button(
                "Close",
                id="close_setting",
                color="dark",
                className="ml-auto",
                outline=True,
            )
        ),
    ],
    id="modal_setting",
    role="document",
    # modalClassName="modal-dialog",
    className="modal-dialog",
)


# end of modal page for advance setting
