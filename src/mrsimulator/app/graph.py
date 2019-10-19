# -*- coding: utf-8 -*-
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc
from mrsimulator.app.toolbar import toolbar

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


spectrum_body = dbc.Card(
    [
        dbc.CardHeader(
            dbc.Row(
                [
                    dbc.Col(
                        html.H4("Spectrum", className="card-title"),
                        # xs=4,
                        # sm=4,
                        # md=4,
                        # lg=4,
                        # xl=4,
                    ),
                    toolbar,
                ]
            )
        ),
        dbc.CardBody([dcc.Graph(id="nmr_spectrum", figure={"data": []})]),
    ],
    className="v-100",
)
