# -*- coding: utf-8 -*-
import dash_core_components as dcc
import dash_html_components as html
import dash_table


__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


colors = {"background": "#e2e2e2", "text": "#585858"}


def top_bar():
    return (
        html.Div(
            className="row",
            children=[
                html.Div(
                    className="col s12 m7 l7",
                    children=[
                        html.H5(id="filename_dataset"),
                        html.H6(
                            id="data_time",
                            style={"textAlign": "left", "color": colors["text"]},
                        ),
                    ],
                ),
                html.Div(
                    className="col s12 m5 l5",
                    children=[
                        dcc.Upload(
                            id="upload_data",
                            children=html.Div(
                                ["Drag and Drop or ", html.A("Select File")]
                            ),
                            style={
                                "width": "100%",
                                "height": "50px",
                                "lineHeight": "50px",
                                "borderWidth": "1px",
                                "borderStyle": "dashed",
                                "borderRadius": "5px",
                                "textAlign": "center",
                                "margin": "1px",
                            },
                            # Allow multiple files to be uploaded
                            multiple=False,
                        ),
                        html.Label(
                            id="error_message",
                            # 'color': 'red'}
                            style={"textAlign": "center"},
                        ),
                    ],
                ),
            ],
        ),
    )


def table_css():
    return {
        "selector": ".dash-cell div.dash-cell-value",
        "rule": {
            "display": "inline",
            "white-space": "inherit",
            "overflow": "inherit",
            "text-overflow": "inherit",
            "overflowX": "scroll",
            "border": "thin lightgrey solid",
        },
    }


def display_isotopomers(isotope, isotopomer_list):
    columns_ = ["index", "site", "isotropic chemical shift", "anisotropy", "asymmetry"]

    child = []
    for j, isotopomer in enumerate(isotopomer_list):
        sites_ = []
        for i, site in enumerate(isotopomer["sites"]):
            if site["isotope"] == isotope:
                site_ = {}
                site_["index"] = i
                site_["isotope"] = site["isotope"]
                site_["isotropic chemical shift"] = site["isotropic_chemical_shift"]
                site_["anisotropy"] = site["shielding_symmetric"]["anisotropy"]
                site_["asymmetry"] = site["shielding_symmetric"]["asymmetry"]

                # ravel_.append(site['isotope'])
                # ravel_.append(site['isotropic_chemical_shift'])
                # ravel_.append(site['shielding_symmetric']['anisotropy'])
                # ravel_.append(site['shielding_symmetric']['asymmetry'])
                sites_.append(site_)

                child.append(
                    html.Div(
                        className="row",
                        children=[
                            html.Div(
                                className="card-panel hoverable col s12 m6 l6",
                                children=[
                                    html.H6(
                                        children="".join(["Isotopomer ", str(j)]),
                                        style={
                                            "textAlign": "left",
                                            "color": colors["text"],
                                        },
                                    ),
                                    dash_table.DataTable(
                                        style_data={"whiteSpace": "normal"},
                                        css=[table_css()],
                                        style_as_list_view=True,
                                        style_cell={
                                            "textAlign": "left",
                                            "padding": "5px",
                                        },
                                        style_header={
                                            "backgroundColor": "white",
                                            "fontWeight": "bold",
                                            "fontSize": 12,
                                            "color": "#585858",
                                        },
                                        style_cell_conditional=[
                                            {
                                                "if": {"row_index": "odd"},
                                                "backgroundColor": "#f8f8f8",
                                            }
                                        ],
                                        columns=[
                                            {"name": i, "id": i} for i in columns_
                                        ],
                                        data=sites_,
                                    ),
                                ],
                            )
                        ],
                    )
                )

    # abundance = isotopomer['abundance']
    return child


def plot_object_widget():
    return [
        html.Div(
            className="card-panel hoverable",
            children=[
                html.H4(
                    id="spectrum_id",
                    children="Spectrum",
                    style={"textAlign": "left", "color": colors["text"]},
                ),
                # dcc.Dropdown(
                #     id="download",
                #     options=[
                #         {"label": "csv", "value": "csv"},
                #         {"label": "csdf", "value": "csdf"},
                #     ],
                #     value="csdf",
                # ),
                html.A(id="download_csv", children="\u21E9 Download CSV"),
                html.A(id="download_csdf", children="\u21E9 Download CSDM"),
                dcc.Graph(id="nmr_spectrum", figure={"data": []}),
            ],
            # style={'width': '100%', 'margin': '0px'}
        )
    ]


def spectrum_object_widget(object_=[]):
    """
    Return the layout for the isotope, number of points,
    spectral width and reference offset.
    """
    return [
        html.Div(
            className="card-panel hoverable",
            children=[
                html.H4(
                    children="Spectrum parameters",
                    style={"textAlign": "left", "color": colors["text"]},
                ),
                # dcc.Tabs(
                #     id="tabs",
                #     value='direct_dimension',
                #     children=[
                #         dcc.Tab(
                #             label='Direct dimension',
                #             value='direct_dimension'
                #         ),
                #         # dcc.Tab(label='Tab Two', value='tab-2-example'),
                #     ]
                # ),
                html.Div(id="tabs_content", children=object_),
            ],
            # style={'width': '100%', 'margin': '0px'}
        )
    ]


def direct_dimension_setup():
    return [
        html.H5(
            children="Direct dimension",
            style={"textAlign": "left", "color": colors["text"]},
        ),
        # environment
        html.Div(
            [html.H6("Environment parameters")],
            style={
                "margin-bottom": "5px",
                "margin-top": "35px",
                "color": colors["text"],
            },
        ),
        # isotope
        html.Label("Isotope"),
        html.Div(
            id="isotope_widget_id",
            children=[dcc.Dropdown(id="isotope_id")],
            style={"margin-bottom": "10px", "margin-top": "0px"},
        ),
        # Magnetic flux density
        html.Label(id="Magnetic_flux_density_output_container"),
        html.Div(
            [
                dcc.Slider(
                    id="magnetic_flux_density", min=0, max=30, step=0.1, value=9.4
                )
            ],
            style={"margin-bottom": "10px", "margin-top": "0px"},
        ),
        # Magic angle spinning
        html.Label(id="spinning_frequency_output_container"),
        html.Div(
            className="row",
            children=[
                dcc.Slider(
                    className="col s6 m6 l6",
                    id="spinning_frequency_in_kHz_coarse",
                    min=0.0,
                    max=110,
                    step=5.0,
                    value=0,
                    marks={0: "0 Hz", 50: "50 kHz", 110: "110 kHz"},
                    # style={'textAlign': 'right'}
                ),
                dcc.Slider(
                    className="col s6 m6 l6",
                    id="spinning_frequency_in_kHz_fine",
                    min=0,
                    max=5,
                    step=0.050,
                    value=2.5,
                    marks={2.5: "+2.5 kHz", 5: "+5 kHz"},
                ),
            ],
            style={"margin-bottom": "10px", "margin-top": "0px"},
        ),
        # dimension
        html.Div(
            className="row",
            children=[
                html.H6(className="col s12 m12 l12", children="Dimension parameters"),
                # daq.BooleanSwitch(
                #     id='ppm_switch',
                #     className='col s6 m6 l6',
                #     label='Show ppm',
                #     labelPosition='bottom',
                #     # size=40,
                #     style={
                #         'margin-bottom': '0px',
                #         'margin-top': '5px',
                #         'color': colors['text']
                #     }
                # )
            ],
            style={
                "margin-bottom": "0px",
                "margin-top": "35px",
                "color": colors["text"],
            },
        ),
        # Number of points
        html.Label(id="number_of_points_output_container"),
        html.Div(
            [
                dcc.Slider(
                    id="number_of_points",
                    min=8,
                    max=16,
                    step=1,
                    value=10,
                    marks={
                        8: "",
                        9: "",
                        10: "",
                        11: "",
                        12: "",
                        13: "",
                        14: "",
                        15: "",
                        16: "",
                    },
                )
            ],
            style={"margin-bottom": "10px", "margin-top": "0px"},
        ),
        # Spectral width
        html.Label(id="frequency_bandwidth_output_container"),
        html.Div(
            className="row",
            children=[
                dcc.Slider(
                    className="col s6 m6 l6",
                    id="frequency_bandwidth_coarse",
                    min=0.0,
                    max=1000,
                    step=50,
                    value=100,
                    marks={0: "0 Hz", 500: "0.5 MHz", 1000: "1 MHz"},
                ),
                dcc.Slider(
                    className="col s6 m6 l6",
                    id="frequency_bandwidth_fine",
                    min=0.0,
                    max=50,
                    step=0.050,
                    value=25,
                    marks={0: "", 25: "+25 kHz", 50: "+50 kHz"},
                ),
            ],
            # style={'margin-bottom': '10px', 'margin-top': '0px'}
        ),
        # Reference offset
        html.Label(id="reference_offset_output_container"),
        html.Div(
            className="row",
            children=[
                dcc.Slider(
                    className="col s6 m6 l6",
                    id="reference_offset_coarse",
                    min=0.0,
                    max=100,
                    step=10,
                    value=0,
                    marks={0: "0 Hz", 50: "50 kHz", 100: "100 kHz"},
                ),
                dcc.Slider(
                    className="col s6 m6 l6",
                    id="reference_offset_fine",
                    min=0,
                    max=10,
                    step=0.050,
                    value=0,
                    marks={0: "", 5: "+5 kHz", 10: "+10 kHz"},
                ),
            ],
            # style={'margin-bottom': '10px', 'margin-top': '0px'}
        ),
    ]
