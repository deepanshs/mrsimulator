# -*- coding: utf-8 -*-
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc


__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


csv_button = dbc.Button(
    "CSV",
    id="download_csv_button",
    outline=True,
    color="primary",
    className="mr-7",
    size="sm",
)
csdm_button = dbc.Button(
    "CSDM",
    id="download_csdm_button",
    outline=True,
    color="primary",
    className="mr-7",
    size="sm",
)

info_button = dbc.Button(
    "Info", id="info", outline=True, color="info", className="mr-1", size="sm"
)

show_individual_button = dbc.Button(
    "Decompose",
    id="decompose",
    outline=True,
    color="info",
    className="mr-1",
    size="sm",
    active=True,
)

button_group = dbc.ButtonGroup(
    [
        dbc.CardLink(
            csv_button, href="", id="download_csv", external_link=True, className="mr-1"
        ),
        dbc.CardLink(
            csdm_button,
            href="",
            id="download_csdm",
            external_link=True,
            className="mr-1",
        ),
    ]
)


modal = html.Div(
    [
        dbc.ButtonGroup([show_individual_button]),
        dbc.Modal(
            [
                dbc.ModalHeader("Header"),
                dbc.ModalBody(
                    "When you run a simulation, it's output will be shown here",
                    id="simulation-output",
                ),
                dbc.ModalFooter(
                    dbc.Button(
                        "Close",
                        id="close",
                        color="dark",
                        className="ml-auto",
                        outline=True,
                    )
                ),
            ],
            id="modal",
        ),
    ]
)

spectrum_body = dbc.Card(
    dbc.CardBody(
        [
            dbc.Row(
                [
                    dbc.Col(
                        html.H4("Spectrum", id="spectrum_id", className="card-title"),
                        xs=4,
                        sm=4,
                        md=4,
                        lg=4,
                        xl=4,
                    ),
                    dbc.Col(modal, align="right"),
                    button_group,
                ]
            ),
            dcc.Graph(id="nmr_spectrum", figure={"data": []}),
        ]
    )
)
