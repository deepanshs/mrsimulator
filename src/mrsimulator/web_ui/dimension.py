# -*- coding: utf-8 -*-
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc
import dash_daq as daq
import numpy as np

# dropdown_menu_items_1 = [
#     dbc.DropdownMenuItem("kHz", id="kHz_1"),
#     dbc.DropdownMenuItem("ppm", id="ppm_1"),
# ]

range_num = [8, 10, 12, 14, 16]
list_of_numbers = {i: f"{2 ** i}" for i in range_num}
field_strength = {2: "200 MHz", 4: "400 MHz", 6: "600 MHz", 8: "0.8 GHz", 10: "1 GHz"}


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

    return html.Div([html.Br(), header, *children])


def coordinate_grid_subgroup(i):
    # number of points
    number_of_points = dbc.FormGroup(
        [
            dbc.Row([
                dbc.Col(
                    dbc.Label("Number of points", color="dark", style={"float": "left"}),
                ),
                dbc.Col(
                    id=f"number_points_label-{i}", style={"float": "right"}
                ),
            ]),
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
            dbc.Input(
                value=25.0,
                min=0.0,
                id=f"spectral_width-{i}",
                pattern="^[-+]?[0-9]*\\.?[0-9]+$",
            ),
            dbc.InputGroupAddon("kHz", addon_type="append"),
        ],
        # size="sm",
    )

    # reference offset
    reference_offset = dbc.InputGroup(
        [
            dbc.InputGroupAddon("Reference offset", addon_type="prepend"),
            dbc.Input(
                value=0.0, id=f"reference_offset-{i}", pattern="^[-+]?[0-9]*\\.?[0-9]+$"
            ),
            dbc.InputGroupAddon("kHz", addon_type="append"),
        ],
        # size="sm",
    )

    broaden_range = {0: "0", 2: "2", 4: "4", 6: "6", 8:"8", 10:"10"}

    line_broadening = html.Div(
        dbc.CardBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Row([
                            dbc.Col(
                                dbc.Label("Line Broadening", style={"float": "left"}),
                            ),
                            dbc.Col(
                                dbc.FormText(
                                    id=f"broadening_sigma_label-{i}", style={"float": "right"})
                                ),
                        ]),
                        dcc.Slider(
                            min=0,
                            max=10,
                            step=1,
                            value=0,
                            marks=broaden_range,
                            id=f"broadening_points-{i}",
                        ),
                        # html.Div(id='slider-output')
                    ]
                )
            ]
        )
    )
    broad_label = dbc.Label(
        size="sm", id=f"test_to_value-{i}", style={"color": "#566573"}
    )

    return [
        number_of_points,
        html.Br(),
        spectral_width,
        reference_offset,
        line_broadening,
        broad_label,
    ]


def environment(i):
    # spectrometer frequency

    spectrometer_frequency = dbc.FormGroup(
        [
            dbc.Row(
                [
                    dbc.Col(
                        dbc.Label(
                            "Spectrometer frequency @1H",
                            color="dark",
                            style={"float": "left"},
                        ),
                        width=9,
                    ),
                    dbc.Col(
                        dbc.FormText(
                            id=f"spectrometer_freq_label-{i}", style={"float": "right"}
                        )
                    ),
                ]
            ),
            dcc.Slider(
                min=1,
                max=12,
                step=0.5,
                value=4,
                marks=field_strength,
                id=f"spectrometer_frequency-{i}",
            ),
        ]
    )

    # rotor frequency
    rotor_frequency = dbc.InputGroup(
        [
            dbc.InputGroupAddon("Rotor frequency", addon_type="prepend"),
            dbc.Input(
                value=0.0,
                id=f"rotor_frequency-{i}",
                pattern="^[-+]?[0-9]*\\.?[0-9]+$",
                min=0.0,
                list=["0", "54.7356", "30", "60", "90"],
            ),
            dbc.InputGroupAddon("kHz", addon_type="append"),
        ],
        # size="sm",
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
            dbc.Input(
                value=54.735,
                id=f"rotor_angle-{i}",
                pattern="^[-+]?[0-9]*\\.?[0-9]+$",
                max=90,
                min=0,
                # step=0.01,
                # type="number",
                # inputMode="numeric",
            ),
            dbc.InputGroupAddon("deg", addon_type="append"),
        ],
        # size="sm",
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
                        dbc.Label("Isotope", color="dark"),
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

    return html.Div(
        [
            isotope_and_filter,
            spectrometer_frequency,
            html.Br(),
            rotor_frequency,
            rotor_angle,
        ],
        style={"padding": 0},
    )


# dimension parameters
def make_dimension(i):

    # dimension parameters
    dimension_contents = dbc.Tab(
        label=f"Index-{i}",
        children=[
            environment(i),
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


# badge = dbc.Badge(
#     " ", pill=True, color="success", className="mr-1", id="indicator_status"
# )


dimension_body = dbc.Card(
    [
        dbc.CardBody(
            [
                dbc.Row(dbc.Col(html.H4("Dimensions", className="card-title"))),
                dbc.Tabs([make_dimension(i) for i in range(2)]),
            ],
            className="w-100",
        )
    ],
    className="h-100 my-card-2",
)
