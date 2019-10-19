# -*- coding: utf-8 -*-
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc
from mrsimulator.app.custom_widgets import custom_slider, custom_input_group

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


# dropdown_menu_items_1 = [
#     dbc.DropdownMenuItem("kHz", id="kHz_1"),
#     dbc.DropdownMenuItem("ppm", id="ppm_1"),
# ]


def sub_group(text, id, children=None, hide=True):
    # header = dbc.Row(
    #     [
    #         dbc.Col(html.Hr()),
    #         dbc.Col(dbc.Label(html.Span(text, id=id)), className="col-auto"),
    #         dbc.Col(html.Hr()),
    #     ]
    # )
    if hide:
        classname = "filter-content collapse hide"
    else:
        classname = "filter-content collapse show"

    return html.Div(
        [
            html.Div(
                html.A(
                    [
                        # html.Br(),
                        html.H6(
                            [text, html.I(className="icon-action fas fa-chevron-down")]
                        )
                    ],
                    href="#",
                    **{
                        "data-toggle": "collapse",
                        "data-target": f"#{id}",
                        "aria-expanded": True,
                    },
                    style={"color": "black"},
                ),
                className="my-subcard",
                style={
                    "padding-top": 4,
                    "padding-bottom": 1,
                    "padding-left": 10,
                    "padding-right": 6,
                },
            ),
            html.Div(
                className=classname,
                id=id,
                children=html.Div([*children, html.P()]),
                # style={"padding-bottom": 10},
            ),
        ],
        # style={"border-bottom": "1px solid rgb(214, 214, 214)"},
    )


def coordinate_grid_subgroup(i):
    # number of points
    range_num = [8, 10, 12, 14, 16]
    list_of_numbers = {i: f"{2 ** i}" for i in range_num}
    number_of_points = custom_slider(
        label="Number of points",
        return_function=lambda x: 2 ** x,
        min=7,
        max=17,
        step=1,
        value=11,
        marks=list_of_numbers,
        id=f"number_of_points-{i}",
    )

    # spectral width
    spectral_width = custom_input_group(
        prepend_label="Spectral width",
        append_label="kHz",
        value=25.0,
        min=0.0,
        id=f"spectral_width-{i}",
    )

    # reference offset
    reference_offset = custom_input_group(
        prepend_label="Reference offset",
        append_label="kHz",
        value=0.0,
        id=f"reference_offset-{i}",
    )

    return [number_of_points, html.Br(), spectral_width, reference_offset]


def post_simulation_widgets(i):
    # line broadening ------------------------------------------------------- #
    broaden_range = {0: "0", 2: "2", 4: "4", 6: "6", 8: "8", 10: "10"}
    line_broadening = custom_slider(
        label="Line Broadening",
        return_function=lambda x: f"{x}  ppm",
        min=0,
        max=10,
        step=1,
        value=0,
        marks=broaden_range,
        id=f"broadening_points-{i}",
    )

    return [line_broadening]


def environment(i):
    # spectrometer frequency
    field_strength = {
        2: "200 MHz",
        4: "400 MHz",
        6: "600 MHz",
        8: "0.8 GHz",
        10: "1 GHz",
    }
    spectrometer_frequency = custom_slider(
        label="Spectrometer frequency @1H",
        return_function=lambda x: f"{int(x*100)} MHz" if x < 10 else f"{x/10} GHz",
        min=1,
        max=12,
        step=0.5,
        value=4,
        marks=field_strength,
        id=f"spectrometer_frequency-{i}",
    )

    # rotor frequency
    rotor_frequency = custom_input_group(
        prepend_label="Rotor frequency",
        append_label="kHz",
        value=0.0,
        id=f"rotor_frequency-{i}",
        min=0.0,
        list=["0", "54.7356", "30", "60", "90"],
    )

    # rotor angle
    rotor_angle = custom_input_group(
        prepend_label="Rotor angle",
        append_label="deg",
        value=54.735,
        id=f"rotor_angle-{i}",
        max=90,
        min=0,
    )

    # filter_spin = [
    #     {"label": "1/2", "value": 0.5},
    #     {"label": "1", "value": 1},
    #     {"label": "3/2", "value": 1.5},
    #     {"label": "5/2", "value": 2.5},
    # ]

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
            sub_group("Environment", f"environment_id-{i}", environment(i), hide=False),
            sub_group(
                "Coordinate grid",
                f"coordinate_grid_id-{i}",
                coordinate_grid_subgroup(i),
                hide=False,
            ),
            # sub_group(
            #     "Post simulation", f"post_simulation_id-{i}",
            #      post_simulation_widgets(i)
            # ),
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
                dbc.Tabs([make_dimension(i) for i in range(1)]),
            ],
            className="w-100",
        )
    ],
    className="h-100 my-card-2",
)
