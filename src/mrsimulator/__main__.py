# -*- coding: utf-8 -*-
import base64
import io
import json
import os

import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import numpy as np
import flask
import plotly.graph_objs as go

from copy import deepcopy
from dash import Dash
from dash.dependencies import Input
from dash.dependencies import Output
from dash.dependencies import State

from mrsimulator import Simulator, Isotopomer, SpectroscopicDimension
from mrsimulator.methods import one_d_spectrum
from mrsimulator.web_ui.widgets import get_isotopomers, main_body
from mrsimulator.web_ui import navbar, sidebar

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


external_scripts = [
    {
        "src": "https://use.fontawesome.com/releases/v5.0.13/js/solid.js",
        "integrity": (
            "sha384-tzzSw1/Vo+0N5UhStP3bvwWPq+uvzCMfrN1fEFe+xBmv1C/AtVX5K0uZtmcHitFZ"
        ),
        "crossorigin": "anonymous",
    },
    {
        "src": "https://use.fontawesome.com/releases/v5.0.13/js/fontawesome.js",
        "integrity": (
            "sha384-6OIrr52G08NpOFSZdxxz1xdNSndlD4vdcf/q2myIUVO0VsqaGHJsB0RaBE01VTOY"
        ),
        "crossorigin": "anonymous",
    },
    {
        "src": "https://code.jquery.com/jquery-3.3.1.slim.min.js",
        "integrity": (
            "sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo"
        ),
        "crossorigin": "anonymous",
    },
    {
        "src": (
            "https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.0/umd/popper.min.js"
        ),
        "integrity": (
            "sha384-cs/chFZiN24E4KMATLdqdvsezGxaGsi4hLGOzlXwp5UZB1LY//20VyM2taTB4QvJ"
        ),
        "crossorigin": "anonymous",
    },
    {
        "src": "https://stackpath.bootstrapcdn.com/bootstrap/4.1.0/js/bootstrap.min.js",
        "integrity": (
            "sha384-uefMccjFJAIv6A+rW+L4AHf99KvxDjWSu1z9VI8SKNVmz4sk7buKt/6v9KI65qnm"
        ),
        "crossorigin": "anonymous",
    },
]


app = Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    external_scripts=external_scripts,
    meta_tags=[{"name": "viewport", "content": "width=device-width"}],
)

server = app.server


FIRST = True
sim = Simulator()
sim.x = np.asarray([-1.0, 1.0])
sim.y = np.asarray([0.0, 0.0])
sim.nt = 64
sim.spectrum_previous = {}
sim.new = False
count = 0
total_sites = 0

app.layout = html.Div(
    [
        dbc.Container(
            [
                navbar.navbar,
                dbc.Row(
                    [
                        dbc.Col(sidebar.sidebar, xs=12, sm=12, md=12, lg=12, xl=3),
                        dbc.Col(
                            [
                                html.Br(),
                                html.Div(main_body),
                                html.Br(),
                                html.Div(id="isotopomer_computed_log"),
                            ],
                            xs=12,
                            sm=12,
                            md=12,
                            lg=12,
                            xl=9,
                        ),
                    ]
                ),
                dbc.Jumbotron(),
            ],
            fluid=True,
        )
    ],
    className="flex-display",
)


# toggle download buttons
@app.callback(
    [
        Output("download_csv_button", "disabled"),
        Output("download_csdm_button", "disabled"),
    ],
    [Input("nmr_spectrum", "figure")],
    [State("filename_dataset", "children")],
)
def toggle_download_buttons(value, filename_dataset):
    """Toggle download buttons---csv, csdm."""
    if filename_dataset in [None, ""]:
        return [True, True]
    else:
        return [False, False]


# update the link to the downloadable serialized file.
@app.callback(
    [Output("download_csv", "href"), Output("download_csdm", "href")],
    [
        Input("nmr_spectrum", "figure")
        # Input("download_csv_button", "n_clicks"),
        # Input("download_csdm_button", "n_clicks"),
    ],
    [
        State("filename_dataset", "children"),
        State("isotope_id-0", "value"),
        State("download_csv_button", "n_clicks_timestamp"),
        State("download_csdm_button", "n_clicks_timestamp"),
        State("decompose", "active"),
        State("spinning_frequency-0", "value"),
        State("number_of_points-0", "value"),
        State("frequency_bandwidth-0", "value"),
        State("reference_offset-0", "value"),
        State("magnetic_flux_density-0", "value"),
    ],
)
def update_link(
    figure,
    filename_dataset,
    isotope_id,
    t1,
    t2,
    decompose_status,
    spinning_frequency,
    number_of_points,
    frequency_bandwidth,
    reference_offset,
    magnetic_flux_density,
):
    """Update the link to the downloadable serialized file."""
    name = os.path.splitext(str(filename_dataset))[0]

    global count
    count += 1

    head = "/dash/urlToDownload?name={0}&isotope={1}".format(name, isotope_id)
    tail = "&status={0}&index={1}_{2}_{3}_{4}_{5}_{6}".format(
        str(decompose_status),
        str(spinning_frequency),
        str(number_of_points),
        str(frequency_bandwidth),
        str(reference_offset),
        str(magnetic_flux_density),
        str(count),
    )
    v1 = head + "&id=csv" + tail
    v2 = head + "&id=csdf" + tail
    return [v1, v2]


# Serialize the computed spectrum and download the serialized file.
@app.server.route("/dash/urlToDownload")
def download():
    """Serialize the computed spectrum and download the serialized file."""
    # creating a dynamic csv or file here using `StringIO`
    filename = flask.request.args.get("name")
    isotope_id = flask.request.args.get("isotope")
    file_type = flask.request.args.get("id")
    index = flask.request.args.get("index")
    is_decomposed = flask.request.args.get("status")

    str_io = io.StringIO()

    global count

    mem = io.BytesIO()
    file_name = "_".join([filename, isotope_id, is_decomposed, index, str(count)])

    if file_type == "csv":
        writer = np.asarray([sim.x, sim.y]).T
        header = (("\n {0} from file {1}.json\nfrequency / kHz, amplitudes\n")).format(
            isotope_id, filename
        )

        # save file as csv
        np.savetxt(str_io, writer, fmt="%f", delimiter=",", header=header)

        mem.write(str_io.getvalue().encode("utf-8"))
        mem.seek(0)
        str_io.close()

        return flask.send_file(
            mem,
            mimetype="text/csv",
            attachment_filename=f"{file_name}.csv",
            as_attachment=True,
        )

    if file_type == "csdf":
        new = sim.as_csdm_object()
        new.save(output_device=str_io)
        mem.write(str_io.getvalue().encode("utf-8"))
        mem.seek(0)
        str_io.close()
        return flask.send_file(
            mem, attachment_filename=f"{file_name}.csdf", as_attachment=True
        )


# Main function. Evaluates the spectrum and update the plot.
@app.callback(
    [Output("nmr_spectrum", "figure"), Output("indicator_status", "color")],
    [
        Input("spinning_frequency-0", "value"),
        Input("rotor_angle-0", "value"),
        Input("number_of_points-0", "value"),
        Input("frequency_bandwidth-0", "value"),
        Input("reference_offset-0", "value"),
        Input("magnetic_flux_density-0", "value"),
        Input("isotope_id-0", "value"),
        Input("decompose", "active"),
    ],
)
def update_data(
    spinning_frequency,
    rotor_angle,
    number_of_points,
    frequency_bandwidth,
    reference_offset,
    magnetic_flux_density,
    isotope_id,
    decompose,
):
    """Evaluate the spectrum and update the plot."""
    color = "danger"  # "#ff2b2b"
    if frequency_bandwidth in [None, 0, "", ".", "-"]:
        return [clear_1D_plot(), color]
    if reference_offset in [None, "", ".", "-"]:
        return [clear_1D_plot(), color]
    if spinning_frequency in [None, "", ".", "-"]:
        return [clear_1D_plot(), color]

    # exit when the following conditions are True
    if number_of_points == 0 or isotope_id in ["", None]:
        return [clear_1D_plot(), color]

    # calculating frequency_bandwidth
    try:
        frequency_bandwidth = float(frequency_bandwidth)
    except ValueError:
        return [clear_1D_plot(), color]

    # calculating spin_frequency
    try:
        spin_frequency = float(eval(str(spinning_frequency)))
    except ValueError:
        return [clear_1D_plot(), color]
    except SyntaxError:
        return [clear_1D_plot(), color]

    # calculating reference_offset
    try:
        reference_offset = float(reference_offset)
    except ValueError:
        return [clear_1D_plot(), color]

    try:
        magnetic_flux_density = float(magnetic_flux_density) / 42.57747892
    except ValueError:
        return [clear_1D_plot(), color]

    # if MAS_switch:
    try:
        rotor_angle_in_degree = float(eval(str(rotor_angle)))  # 54.735
    except ValueError:
        return [clear_1D_plot(), color]
    except SyntaxError:
        return [clear_1D_plot(), color]

    spectrum = {
        "isotope": isotope_id,
        "magnetic_flux_density": str(magnetic_flux_density * 100) + " T",
        "rotor_frequency": str(spin_frequency) + " kHz",
        "rotor_angle": str(rotor_angle_in_degree) + " deg",
        "number_of_points": 2 ** number_of_points,
        "spectral_width": str(frequency_bandwidth) + " kHz",
        "reference_offset": str(reference_offset) + " kHz",
    }

    if not sim.spectrum_previous == spectrum or sim.new:
        sim.new = False
        sim.spectrum_previous = deepcopy(spectrum)

        sim.spectrum = [SpectroscopicDimension.parse_dict_with_units(spectrum)]
        sim.x, sim.y = sim.run(
            one_d_spectrum,
            geodesic_polyhedron_frequency=sim.nt,
            individual_spectrum=True,
        )

        sim.x = sim.x.value

    return plot_1D(isotope_id, decompose)


# Update the isotopomers when a new file is imported.
@app.callback(
    [
        Output("filename_dataset", "children"),
        Output("data_time", "children"),
        Output("error_message", "children"),
        Output("isotope_id-0", "options"),
        Output("isotope_id-0", "value"),
    ],
    [Input("upload_data", "contents")],
    [State("upload_data", "filename"), State("upload_data", "last_modified")],
)
def update_isotopomers(content, filename, date):
    """Update the isotopomers when a new file is imported."""
    children, success = parse_contents(content, filename, date)

    if success:
        isotope_list = [
            {"label": site_iso, "value": site_iso} for site_iso in sim.get_isotopes()
        ]
        isotope = isotope_list[0]["value"]
        sim.new = True
    else:
        sim.isotopomers = []
        isotope_list = []
        isotope = None
    return [children[0], children[1], children[2], isotope_list, isotope]


# Model window for advance input,
@app.callback(
    Output("modal_setting", "is_open"),
    [Input("advance-id", "n_clicks"), Input("close_setting", "n_clicks")],
    [State("modal_setting", "is_open")],
)
def toggle_modal_setting(n1, n2, is_open):
    """Model window for advance input."""
    if n1 or n2:
        return not is_open
    return is_open


# Model window option, number-of-orientations, from advance input.
@app.callback(
    Output("number_of_averaging_points", "children"),
    [Input("averaging_quality", "value")],
)
def update_number_of_orientations(value):
    """
    Update the number of orientation for powder averaging.
    Option for advance modal.
    """
    sim.nt = value
    ori = 2 * (value + 1) * (value + 2)
    return f"Averaging over {ori} orientations.".format(value)


# @app.callback(
#     Output("unit_sw", "label"),
#     [
#       Input("Hz_1", "n_clicks"),
#       Input("kHz_1", "n_clicks"),
#       Input("MHz_1", "n_clicks")
#     ],
# )
# def update_label_1(n1, n2, n3):
#     id_lookup = {"Hz_1": "Hz", "kHz_1": "kHz", "MHz_1": "MHz"}
#     return update_label(n1, n2, n3, id_lookup)


# @app.callback(
#     Output("unit_ref", "label"),
#     [
#       Input("Hz_2", "n_clicks"),
#       Input("kHz_2", "n_clicks"),
#       Input("MHz_2", "n_clicks")
#     ],
# )
# def update_label_2(n1, n2, n3):
#     id_lookup = {"Hz_2": "Hz", "kHz_2": "kHz", "MHz_2": "MHz"}
#     return update_label(n1, n2, n3, id_lookup)


# @app.callback(
#     Output("unit_spin", "label"),
#     [
#       Input("Hz_3", "n_clicks"),
#       Input("kHz_3", "n_clicks"),
#       Input("MHz_3", "n_clicks")
#     ],
# )
# def update_label_3(n1, n2, n3):
#     id_lookup = {"Hz_3": "Hz", "kHz_3": "kHz", "MHz_3": "MHz"}
#     return update_label(n1, n2, n3, id_lookup)


# def update_label(n1, n2, n3, id_lookup):
#     # use a dictionary to map ids back to the desired label
#     # makes more sense when there are lots of possible labels
#     # id_lookup = {"Hz": "Hz", "kHz": "kHz", "MHz": "MHz"}

#     ctx = dash.callback_context

#     if (n1 is None and n2 is None and n3 is None) or not ctx.triggered:
#         # if neither button has been clicked, return "Not selected"
#         return "kHz"

#     # this gets the id of the button that triggered the callback
#     button_id = ctx.triggered[0]["prop_id"].split(".")[0]
#     return id_lookup[button_id]


# add callback for toggling the collapse on small screens
# @app.callback(
#     Output("navbar-collapse", "is_open"),
#     [Input("navbar-toggler", "n_clicks")],
#     [State("navbar-collapse", "is_open")],
# )
# def toggle_navbar_collapse(n, is_open):
#     if n:
#         return not is_open
#     return is_open


# @app.callback(
#     [
#         Output(f"collapse-dimension-0", "is_open"),
#         # *[Output(f"collapse-card-{i}", "n_clicks") for i in range(total_sites)],
#     ],
#     [
#         Input(f"dimension-toggle-0", "n_clicks"),
#         # *[Input(f"card-toggle-{i}", "n_clicks") for i in range(total_sites)],
#     ],
#     [
#         State(f"collapse-dimension-0", "is_open"),
#         # *[State(f"collapse-card-{i}", "n_clicks") for i in range(total_sites)],
#     ],
# )
# def toggle_accordion(n1, is_dimension_open):
#     ctx = dash.callback_context
#     # global total_sites
#     # print("total_sites", total_sites)
#     # list_output = [False for _ in range(total_sites + 1)]

#     if not ctx.triggered:
#         return [""]
#     else:
#         button_id = ctx.triggered[0]["prop_id"].split(".")[0]

#     if button_id == "dimension-toggle-0" and n1:
#         # list_output[0] = not is_dimension_open
#         return [not is_dimension_open]
#     # for i in range(total_sites):
#     #     if button_id == "card-toggle-{i}" and n_list[i]:
#     #         list_output[i] = not is_open[i]
#     #     return False, not is_open2, False
#     # elif button_id == "group-3-toggle" and n3:
#     #     return False, False, not is_open3
#     return [False]


@app.callback(
    [Output("decompose", "active")],
    [Input("decompose", "n_clicks")],
    [State("decompose", "active")],
)
def toggle_decompose_button(n1, status):
    "Toggle decompose button."
    new_status = True
    if bool(status):
        new_status = False
    return [new_status]


# @app.callback(
#     [Output("isotopomer_computed_log", "children")],
#     [Input("isotope_id", "value"), Input("magnetic_flux_density", "value")],
# )
# def update_isotopomer_list(isotope, flux):
#     cards_deck = get_isotopomers(isotope, sim.isotopomers)

#     return [cards_deck]


# @app.callback(
#     [Output("dimension-body", "children")], [Input("dimension-tabs", "active_tab")]
# )
# def tab_content(active_tab):
#     index = int(active_tab.split("-")[1])
#     return [dimension_children[index]]


def clear_1D_plot():
    """Create and return a new blank plot."""
    data = go.Scatter(
        x=[-1, 1],
        y=[0, 0],
        text="",
        mode="lines",
        opacity=1.0,
        line={"color": "black", "width": 1.2},
    )
    data_object = {
        "data": [data],
        "layout": go.Layout(
            xaxis={
                "type": "linear",
                "title": "frequency / kHz",
                "ticks": "outside",
                "showline": True,
                "autorange": True,
            },
            yaxis={
                "type": "linear",
                "title": "arbitrary unit",
                "ticks": "outside",
                "showline": True,
                "autorange": True,
            },
            autosize=True,
            transition={"duration": 500},
            margin={"l": 50, "b": 40, "t": 5, "r": 5},
            legend={"x": 0, "y": 1},
            hovermode="closest",
        ),
    }
    return data_object


def plot_1D(isotope_id, decompose):
    """Generate and return a one-dimensional plot instance."""
    data = []
    if decompose:
        for i, datum in enumerate(sim.y):
            if not isinstance(datum, list):
                name = sim.isotopomers[i].name
                if name == "":
                    name = "Isotopomer " + str(i + 1)
                data.append(
                    go.Scatter(
                        x=sim.x,
                        y=datum,
                        mode="lines",
                        opacity=0.8,
                        line={"width": 1.2},
                        fill="tozeroy",
                        name=name,
                        # hovername=
                        # customdata=str(sim.isotopomers[i]),
                    )
                )
    else:
        y_data = 0
        for i, datum in enumerate(sim.y):
            if not isinstance(datum, list):
                y_data += datum
        data.append(
            go.Scatter(
                x=sim.x,
                y=y_data,
                mode="lines",
                line={"color": "black", "width": 1.2},
                opacity=1.0,
                # name=isotope_id,
                name=f"spectrum",
                # hovername=
                # customdata=str(sim.isotopomers[i]),
            )
        )

    x_label = str(isotope_id + f" frequency / ppm")

    data_object = {
        "data": data,
        "layout": go.Layout(
            xaxis=dict(
                type="linear",
                title=x_label,
                ticks="outside",
                showline=True,
                autorange="reversed",
                zeroline=False,
            ),
            yaxis=dict(
                type="linear",
                title="arbitrary unit",
                ticks="outside",
                showline=True,
                zeroline=False,
                rangemode="tozero",
            ),
            autosize=True,
            transition={"duration": 150, "easing": "sin-out"},
            margin={"l": 50, "b": 40, "t": 5, "r": 5},
            legend={"x": 0, "y": 1},
            hovermode="closest",
        ),
    }
    return [data_object, "success"]


# line={"shape": "hv", "width": 1},

#  ['linear', 'quad', 'cubic', 'sin', 'exp', 'circle',
#             'elastic', 'back', 'bounce', 'linear-in', 'quad-in',
#             'cubic-in', 'sin-in', 'exp-in', 'circle-in', 'elastic-in',
#             'back-in', 'bounce-in', 'linear-out', 'quad-out',
#             'cubic-out', 'sin-out', 'exp-out', 'circle-out',
#             'elastic-out', 'back-out', 'bounce-out', 'linear-in-out',
#             'quad-in-out', 'cubic-in-out', 'sin-in-out', 'exp-in-out',
#             'circle-in-out', 'elastic-in-out', 'back-in-out',
#             'bounce-in-out']


def parse_contents(contents, filename, date):
    """Parse contents from the isotopomers file."""
    try:
        if "json" in filename:
            content_string = contents.split(",")[1]
            decoded = base64.b64decode(content_string)
            data = json.loads(str(decoded, encoding="UTF-8"))
            isotopomers = data["isotopomers"]

            name = filename
            if "name" in data.keys():
                name = data["name"]
                if name == "":
                    name = filename

            description = ""
            if "description" in data.keys():
                description = data["description"]

            sim.isotopomers = [
                Isotopomer.parse_dict_with_units(item) for item in isotopomers
            ]
            sim.spectrum = []
            return (
                [name, description, "Select a JSON serialized .isotopomers file."],
                True,
            )

        else:
            return (
                ["", "", "A JSON file with valid list of isotopomers is required."],
                False,
            )

    except Exception:
        if FIRST:
            return (["", "", "Select a JSON serialized .isotopomers file."], False)
            # FIRST = False
        else:
            return ["", "", "Error reading the file."], False


if __name__ == "__main__":
    app.run_server(debug=True)
