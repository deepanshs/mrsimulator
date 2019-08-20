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
    color="info",
    className="mr-1",
    size="sm",
)
csdm_button = dbc.Button(
    "CSDM",
    id="download_csdm_button",
    outline=True,
    color="info",
    className="mr-1",
    size="sm",
)

info_button = dbc.Button(
    "Info", id="info", outline=True, color="info", className="mr-1", size="sm"
)

button_group = dbc.ButtonGroup(
    [
        dbc.CardLink(csv_button, href="", id="download_csv", external_link=True),
        dbc.CardLink(csdm_button, href="", id="download_csdm", external_link=True),
    ]
)


modal = html.Div(
    [
        dbc.ButtonGroup([info_button]),
        dbc.Modal(
            [
                dbc.ModalHeader("Header"),
                dbc.ModalBody(
                    "When you run a simulation, it's output will be shown here",
                    id="simulation-output",
                ),
                dbc.ModalFooter(
                    dbc.Button("Close", id="close", className="ml-auto", outline=True)
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
                        html.H3("Spectrum", id="spectrum_id", className="card-title"),
                        xs=4,
                        sm=4,
                        md=4,
                        lg=4,
                        xl=4,
                    ),
                    dbc.Col(modal, align="right"),
                    button_group,
                    #     ]
                    # ),
                    # # html.H3("Spectrum", id="spectrum_id", className="card-title"),
                    # dbc.Col(
                    #     [
                    # dbc.Col(
                    # dbc.Button(
                    #     dbc.DropdownMenuItem(
                    #         "CSV",
                    #         id="download_csv",
                    #         external_link=True,
                    #         style={
                    #             "background-color": "rgba(245, 245, 245, 0)",
                    #             "opacity": "1.0",
                    #             "border": "none",
                    #         },
                    #     ),
                    #     outline=True,
                    #     color="info",
                    #     className="mr-1",
                    #     size="sm",
                    # ),
                    # ),
                    # dbc.Col(
                    # dbc.Button(
                    #     dbc.DropdownMenuItem(
                    #         "CSDM", id="download_csdf", external_link=True
                    #     ),
                    #     # id="download_csdm",
                    #     outline=True,
                    #     color="info",
                    #     className="mr-1",
                    #     size="sm",
                    # ),
                    # ),
                    # dbc.Col(html.H5("info", id="spectrum_id_info")),
                    # dbc.DropdownMenu(
                    #     [
                    #         # dbc.DropdownMenuItem(
                    #         #     "CSV", id="download_csv", external_link=True
                    #         # ),
                    #         # dbc.DropdownMenuItem(
                    #         #     "CSDM", id="download_csdf", external_link=True
                    #         # )
                    #     ],
                    #     label="Download",
                    #     group=True,
                    # ),
                    # html.A("\u21E9 Download CSV", id="download_csv"),
                    # html.A("\u21E9 Download CSDM", id="download_csdf"),
                ]
            ),
            # html.P(id="placeholder"),
            html.Div(html.A(id="download")),
            dcc.Graph(id="nmr_spectrum", figure={"data": []}),
        ]
    )
)
