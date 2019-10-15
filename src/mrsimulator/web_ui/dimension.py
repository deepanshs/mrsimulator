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
field_strength = {2: "200 MHz", 4: "400 MHz", 6: "600 MHz", 8: "800 MHz", 10: "1 GHz"}


def sub_group(text, id, children):
    header = dbc.Row(
        [
            dbc.Col(html.Hr()),
            dbc.Col(dbc.Label(html.Span(text, id=id)), className="col-auto"),
            dbc.Col(html.Hr()),
        ]
    )

    # return html.Div(
    #     className="card card-info",
    #     children=[
    #         html.Div(
    #             className="my-card-header",
    #             children=[
    #                 html.Span(
    #                     html.H6(
    #                         [text, html.I(className="fas fa-chevron-down")],
    #                         className="card-title",
    #                     ),
    #                     id=id,
    #                     # className="pull-right clickable",
    #                 )
    #             ],
    #         ),
    #         html.Div(className="card-body", children=children),
    #     ],
    # )

    return html.Div([header, *children])


def coordinate_grid_subgroup(i):
    # number of points
    number_of_points = dbc.FormGroup(
        [
            dbc.Label("Number of points"),
            dcc.Slider(
                min=7,
                max=17,
                step=1,
                value=11,
                marks=list_of_numbers,
                id=f"number_of_points-{i}",
            ),
        ]
    )

    # spectral width
    spectral_width = dbc.InputGroup(
        [
            dbc.InputGroupAddon("Spectral width", addon_type="prepend"),
            dbc.Input(value=25.0, min=0.0, id=f"frequency_bandwidth-{i}"),
            dbc.InputGroupAddon("kHz", addon_type="append"),
        ]
    )

    # reference offset
    reference_offset = dbc.InputGroup(
        [
            dbc.InputGroupAddon("Reference offset", addon_type="prepend"),
            dbc.Input(value=0.0, id=f"reference_offset-{i}"),
            dbc.InputGroupAddon("kHz", addon_type="append"),
        ]
    )

    return [number_of_points, html.Br(), spectral_width, reference_offset]


def environment(i):
    # spectrometer frequency

    spectrometer_frequency = dbc.FormGroup(
        [
            dbc.Label("Spectrometer frequency @ 1H"),
            dcc.Slider(
                min=1,
                max=11,
                step=0.5,
                value=4,
                marks=field_strength,
                id=f"magnetic_flux_density-{i}",
            ),
        ]
    )

    # rotor frequency
    rotor_frequency = dbc.InputGroup(
        [
            dbc.InputGroupAddon("Rotor frequency", addon_type="prepend"),
            dbc.Input(value=0, id=f"spinning_frequency-{i}"),
            dbc.InputGroupAddon("kHz", addon_type="append"),
        ]
    )

    # rotor angle
    # rotor_marks = {
    #     0: "0°",
    #     10: "10°",
    #     20: "20°",
    #     30: "30°",
    #     40: "40°",
    #     50: "50°",
    #     54.735: "MA",
    #     60: "60°",
    #     70: "70°",
    #     80: "80°",
    #     90: "90°",
    # }
    # rotor_angle = dbc.FormGroup(
    #     [
    #         dbc.Label("Rotor angle"),
    #         dcc.Slider(
    #             min=0, max=90, step=10, value=54.735,
    #             marks=rotor_marks, id="rotor_angle"
    #         ),
    #     ]
    # )

    rotor_angle = dbc.InputGroup(
        [
            dbc.InputGroupAddon("Rotor angle", addon_type="prepend"),
            dbc.Input(value=54.735, id=f"rotor_angle-{i}"),
            dbc.InputGroupAddon("degree", addon_type="append"),
        ]
    )

    filter_spin = [
        {"label": "1/2", "value": 0.5},
        {"label": "1", "value": 1},
        {"label": "3/2", "value": 1.5},
        {"label": "5/2", "value": 2.5},
    ]

    isotope_and_filter = dbc.Row(
        [
            dbc.Col(
                dbc.FormGroup(
                    [
                        dbc.Label("Isotope", className="mr-2"),
                        dcc.Dropdown(id=f"isotope_id-{i}"),
                    ]
                )
            ),
            # dbc.Col(
            #     dbc.FormGroup(
            #         [
            #             dbc.Label("Filter", className="mr-2"),
            #             dcc.Dropdown(id="filter_spin", options=filter_spin,
            #                          value=0.5),
            #         ]
            #     )
            # ),
        ]
    )

    return [
        isotope_and_filter,
        spectrometer_frequency,
        html.Br(),
        rotor_frequency,
        rotor_angle,
    ]


# dimension parameters
def make_dimension(i):

    # dimension parameters
    dimension_contents = dbc.Tab(
        label=f"Index-{i}",
        children=[
            *environment(i),
            # html.Br(),
            sub_group(
                "Coordinate grid",
                f"coordinate_grid_id={i}",
                coordinate_grid_subgroup(i),
            ),
        ],
        # style={
        #     "min-height": "65vh",
        #     "max-height": "65vh",
        #     "overflow-y": "scroll",
        #     "overflow-x": "hidden",
        # },
    )
    return dimension_contents


# submit_button = dbc.Button(
#     "Submit",
#     id="submit_query",
#     outline=True,
#     color="primary",
#     className="mr-1",
#     size="sm",
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
    [
        dbc.CardBody(
            [
                dbc.Row(
                    [
                        dbc.Col(html.H4("Dimensions", className="card-title")),
                        dbc.Col(modal),
                    ]
                ),
                dbc.Tabs([make_dimension(i) for i in range(2)]),
            ],
            className="w-100",
        )
    ],
    className="h-100 my-card-2",
)
