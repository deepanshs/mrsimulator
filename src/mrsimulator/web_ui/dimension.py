# -*- coding: utf-8 -*-
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc
import dash_daq as daq
import numpy as np

dropdown_menu_items_1 = [
    dbc.DropdownMenuItem("kHz", id="kHz_1"),
    dbc.DropdownMenuItem("ppm", id="ppm_1"),
]

# dimension parameters
range_num = [8, 10, 12, 14, 16]
list_of_numbers = {i: f"{2 ** i}" for i in range_num}

dimension = [
    dbc.FormGroup(
        [
            dbc.Label("Number of points"),
            dcc.Slider(
                min=7,
                max=17,
                step=1,
                value=11,
                marks=list_of_numbers,
                id="number_of_points",
            ),
        ]
    ),
    html.Br(),
    dbc.InputGroup(
        [
            dbc.InputGroupAddon("Spectral width", addon_type="prepend"),
            dbc.Input(value=25.0, min=0.0, id="frequency_bandwidth"),
            dbc.InputGroupAddon("kHz", addon_type="append"),
        ]
    ),
    # html.Br(),
    dbc.InputGroup(
        [
            dbc.InputGroupAddon("Reference offset", addon_type="prepend"),
            dbc.Input(value=0.0, id="reference_offset"),
            dbc.InputGroupAddon("kHz", addon_type="append"),
        ]
    ),
]

collapsible_dimension = dbc.Card(
    [
        dbc.CardHeader(
            [
                # dbc.Label("Dimension parameters"),
                dbc.Button(
                    "Dimension parameters",
                    color="link",
                    size="sm",
                    id="dimension-toggle",
                )
            ]
        ),
        dbc.Collapse(dbc.CardBody(dimension), id="collapse-dimension"),
    ],
    # style={"border": "1px"},
)


# spin and environment parameters

indexes = [2, 4, 6, 8, 10]
field_strength = {2: "200 MHz", 4: "400 MHz", 6: "600 MHz", 8: "800 MHz", 10: "1 GHz"}

filter_spin = [
    {"label": "1/2", "value": 0.5},
    {"label": "1", "value": 1},
    {"label": "3/2", "value": 1.5},
    {"label": "5/2", "value": 2.5},
]

isotope_and_filter = [
    dbc.Col(
        dbc.FormGroup(
            [dbc.Label("Isotope", className="mr-2"), dcc.Dropdown(id="isotope_id")]
        )
    ),
    dbc.Col(
        dbc.FormGroup(
            [
                dbc.Label("Filter", className="mr-2"),
                dcc.Dropdown(id="filter_spin", options=filter_spin, value=0.5),
            ]
        )
    ),
]

spectrometer_frequency = dbc.FormGroup(
    [
        dbc.Label("Spectrometer frequency @ 1H"),
        dcc.Slider(
            min=1,
            max=11,
            step=0.5,
            value=4,
            marks=field_strength,
            id="magnetic_flux_density",
        ),
    ]
)

environment = [
    dbc.Row(isotope_and_filter),
    spectrometer_frequency,
    html.Br(),
    dbc.InputGroup(
        [
            dbc.InputGroupAddon("Spinning frequency", addon_type="prepend"),
            dbc.Input(value=0, id="spinning_frequency"),
            dbc.InputGroupAddon("kHz", addon_type="append"),
        ]
    ),
]


# submit_button = dbc.Button(
#     "Submit",
#     id="submit_query",
#     outline=True,
#     color="primary",
#     className="mr-1",
#     size="sm",
# )

# dimension_body = dbc.Card(
#     dbc.CardBody(
#         [
#             html.H3("Spectrum parameters"),
#             html.H4("Direct dimension"),
#             html.Br(),
#             html.Br(),
#             html.H5("Environment parameters"),
#             html.Div(environment),
#             html.Br(),
#             html.Br(),
#             html.H5("Dimension parameters"),
#             html.Div(dimension),
#             html.Br(),
#             # html.Div(submit_button),
#         ]
#     ),
#     className="mb-0",
# )


badge = dbc.Badge(
    " ", pill=True, color="success", className="mr-1", id="indicator_status"
)

setting_button = dbc.Button(
    "Advance", id="advance-id", color="dark", outline=True, className="mr-1", size="sm"
)

integration_level = {4: "4", 5: "5", 6: "6", 7: "7", 8: "8", 9: "9", 10: "10", 11: "11"}

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

model_line_2 = dbc.Row(
    [
        dbc.Col(dbc.Label("Integration volume")),
        dbc.Col(
            dcc.Dropdown(
                id="octants",
                options=[
                    {"label": "Octant", "value": 0},
                    {"label": "Hemisphere", "value": 1},
                    {"label": "Sphere", "value": 2},
                ],
                value=0,
                clearable=False,
            )
        ),
    ]
)

model_info = dbc.Label(
    size="sm", id="number_of_averaging_points", style={"color": "#566573"}
)

modal = html.Div(
    [
        dbc.ButtonGroup([setting_button, badge]),
        dbc.Modal(
            [
                dbc.ModalHeader("Advance setting"),
                dbc.ModalBody(dbc.FormGroup([model_line_1, model_info])),
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
        ),
    ]
)


dimension_body = dbc.Card(
    dbc.CardBody(
        [
            dbc.Row([dbc.Col(html.H4("Parameters")), dbc.Col(modal, align="right")]),
            html.Div(environment),
            html.Br(),
            html.Div(collapsible_dimension),
        ]
    )
)
