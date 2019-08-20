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

range_num = [8, 10, 12, 14, 16]
list_of_numbers = {i: f"{2 ** i}" for i in range_num}

# data_list = html.Datalist(id="list_of_numbers", children=list_of_numbers)


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
            )
            # dcc.Dropdown(
            #     value=2048,
            #     options=list_of_numbers,
            #     clearable=False,
            #     # min=0,
            #     # max=1000000,
            #     # type="number",
            #     id="number_of_points",
            # )
            # dbc.InputGroupAddon(
            #     dcc.Dropdown(
            #         options=list_of_numbers,
            #         value=2048,
            #         id="number_of_points",
            #         clearable=False,
            #         # style={"width": "20em", "height": "2.37em"},
            #     ),
            #     addon_type="append",
            # ),
            # dbc.Input(
            #     value=2048, min=0, max=1000000, type="number", id="number_of_points"
            # ),
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
    dbc.InputGroup(
        [
            dbc.InputGroupAddon("Reference offset", addon_type="prepend"),
            dbc.Input(value=0.0, id="reference_offset"),
            dbc.InputGroupAddon("kHz", addon_type="append"),
            # dbc.DropdownMenu(
            #     dropdown_menu_items_1,
            #     label="kHz",
            #     addon_type="append",
            #     id="unit_ref",
            # ),
        ]
    ),
]


indexes = [2, 4, 6, 8, 10]
field_strength = {2: "200 MHz", 4: "400 MHz", 6: "600 MHz", 8: "800 MHz", 10: "1 GHz"}

environment = [
    dbc.FormGroup(
        [
            dbc.Label("Isotope", className="mr-2"),
            dcc.Dropdown(
                id="isotope_id",
                # clearable=False,
                # style={"width": "8em", "height": "2.37em"},
            ),
            # dbc.DropdownMenu(
            #     label="",
            #     addon_type="prepend",
            #     id="isotope_id",
            # ),
            # dcc.Dropdown(id="isotope_id"),
        ]
    ),
    # html.Br(),
    dbc.FormGroup(
        [
            dbc.Label("Spectrometer frequency @ 1H"),
            dcc.Slider(
                min=1,
                max=11,
                step=0.5,
                value=4,
                marks=field_strength,
                id="magnetic_flux_density",
            )
            # dcc.Dropdown(
            #     value=400,
            #     options=field_strength,
            #     id="magnetic_flux_density",
            #     clearable=False,
            #     # style={"width": "8em", "height": "2.37em"},
            # ),
            # dbc.Input(value=9.4, min=0, max=100, id="magnetic_flux_density"),
            # dbc.Label("T"),
        ]
    ),
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

indicator = daq.Indicator(id="indicator_status", color="#00cc96", value=True)

dimension_body = [
    dbc.Row(
        [
            dbc.Col(html.H3("Spectrum parameters"), xs=10, sm=10, md=10, lg=10, xl=10),
            dbc.Col(indicator),
        ]
    ),
    html.H4("Direct dimension"),
    # html.Br(),
    # html.Br(),
    html.H5("Environment parameters"),
    html.Div(environment),
    html.Br(),
    # html.Br(),
    html.H5("Dimension parameters"),
    html.Div(dimension),
    # html.Br(),
    # html.Div(submit_button),
]
