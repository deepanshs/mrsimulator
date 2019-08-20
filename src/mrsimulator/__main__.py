# -*- coding: utf-8 -*-
import base64
import datetime
import io
import json
import os

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import numpy as np
import flask
import plotly.graph_objs as go
from dash import Dash
from dash.dependencies import Input
from dash.dependencies import Output
from dash.dependencies import State

from mrsimulator import Simulator, Isotopomer, Spectrum
from mrsimulator.methods import one_d_spectrum
from mrsimulator.web_ui.widgets import display_isotopomers, input_file, main_body

from mrsimulator.unit import string_to_quantity, is_physical_quantity

import csdmpy as cp

from mrsimulator.web_ui import navbar, titlebar


# class Data:
#     def __init__(self):
#         self.x = np.asarray([-1.0, 1.0])
#         self.y = np.asarray([0.0, 0.0])


__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


external_stylesheets = [
    "https://cdnjs.cloudflare.com/ajax/libs/materialize/1.0.0/css/materialize.min.css"
]

external_scripts = [
    "https://cdnjs.cloudflare.com/ajax/libs/materialize/1.0.0/js/materialize.min.js"
]


app = Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    external_scripts=external_scripts,
    # external_stylesheets=external_stylesheets,
)


FIRST = True
sim = Simulator()
sim.x = np.asarray([-1.0, 1.0])
sim.y = np.asarray([0.0, 0.0])
count = 0

app.layout = html.Div(
    className="container",
    style={
        "color": "#585858",
        "background-image": "linear-gradient(#fdfefe, #f4f6f7, #ffffff)",
        # "width": "90%",
    },
    children=[
        navbar.navbar,
        titlebar.titlebar,
        # dcc.ConfirmDialog(
        #     id="confirm", message="Cannot process the current request."
        # ),
        *input_file,
        html.Hr(className="my-2"),
        html.Div(main_body),
        html.Hr(className="my-2"),
        html.Div(id="isotopomer_computed_log"),
    ],
)


# @app.callback([Output("download", "href")], [Input("download_csv", "n_clicks")])
# def update_link(n1):
#     print(n1)


@app.callback(
    [
        Output("download_csv_button", "disabled"),
        Output("download_csdm_button", "disabled"),
    ],
    [Input("nmr_spectrum", "figure")],
    [State("filename_dataset", "children"), State("spectrum_id", "children")],
)
def toggle_download(value, filename_dataset, spectrum_id):
    if filename_dataset in [None, ""] or spectrum_id in [None, ""]:
        return [True, True]
    else:
        return [False, False]


@app.callback(
    [Output("download_csv", "href"), Output("download_csdm", "href")],
    [
        Input("nmr_spectrum", "figure")
        # Input("download_csv_button", "n_clicks"),
        # Input("download_csdm_button", "n_clicks"),
    ],
    [
        State("filename_dataset", "children"),
        State("isotope_id", "value"),
        State("download_csv_button", "n_clicks_timestamp"),
        State("download_csdm_button", "n_clicks_timestamp"),
    ],
)
def update_link(figure, filename_dataset, isotope_id, t1, t2):
    """Update the csv download link when the plot is refreshed."""
    name = os.path.splitext(str(filename_dataset))[0]

    global count
    count += 1
    # print("t1", t1, "t2", t2, "count", count, "isotope_id", isotope_id)

    v1 = "/dash/urlToDownload?name={0}&isotope={1}&id=csv&index={2}".format(
        name, isotope_id, count
    )
    v2 = "/dash/urlToDownload?name={0}&isotope={1}&id=csdf&index={2}".format(
        name, isotope_id, count
    )
    # print(v1, v2)
    return [v1, v2]


@app.server.route("/dash/urlToDownload")
def download():
    # creating a dynamic csv or file here using `StringIO`
    filename = flask.request.args.get("name")
    isotope_id = flask.request.args.get("isotope")

    file_type = flask.request.args.get("id")
    index = flask.request.args.get("index")

    # print(filename, isotope_id, file_type, index)
    str_io = io.StringIO()

    global count

    if file_type == "csv":
        writer = np.asarray([sim.x, sim.y]).T
        _header_ = (
            ("\n {0} from file {1}.json\nfrequency / kHz, amplitudes\n")
        ).format(isotope_id, filename)

        # save file as csv
        np.savetxt(str_io, writer, fmt="%f", delimiter=",", header=_header_)

        mem = io.BytesIO()
        mem.write(str_io.getvalue().encode("utf-8"))
        mem.seek(0)
        str_io.close()
        file_name = "_".join([filename, isotope_id])
        return flask.send_file(
            mem,
            mimetype="text/csv",
            attachment_filename=f"{file_name}.csv",
            as_attachment=True,
        )
    # print(sim.spectrum)
    if file_type == "csdf":
        new = cp.new()
        dimension = {
            "type": "linear",
            "count": sim.spectrum.number_of_points,
            "increment": "{0} Hz".format(
                sim.spectrum.spectral_width / sim.spectrum.number_of_points
            ),
            "coordinates_offset": "{0} Hz".format(sim.spectrum.reference_offset),
            "complex_fft": True,
        }
        dependent_variable = {
            "type": "internal",
            "quantity_type": "scalar",
            "numeric_type": "float64",
            "components": [sim.y],
        }
        new.add_dimension(dimension)
        new.add_dependent_variable(dependent_variable)
        new.dependent_variables[0].encoding = "base64"

        new.save(output_device=str_io)

        mem = io.BytesIO()
        mem.write(str_io.getvalue().encode("utf-8"))
        mem.seek(0)
        str_io.close()
        file_name = "_".join([filename, isotope_id])
        return flask.send_file(
            mem, attachment_filename=f"{file_name}.csdf", as_attachment=True
        )


@app.callback(
    Output("isotopomer_computed_log", "children"), [Input("isotope_id", "value")]
)
def update_isotopomers_table_log(value):
    return display_isotopomers(value, sim.isotopomers)


# @app.callback(
#     Output('tabs_content', 'children'),
#     [Input('tabs', 'value')]
# )
# def update_tab(value):
#     if value == 'direct_dimension':
#         return direct_dimension_setup()


@app.callback(
    [Output("nmr_spectrum", "figure"), Output("indicator_status", "color")],
    # [Input("submit_query", "n_clicks")],
    [
        # Input('confirm', 'submit_n_clicks'),
        Input("spinning_frequency", "value"),
        Input("number_of_points", "value"),
        Input("frequency_bandwidth", "value"),
        Input("reference_offset", "value"),
        Input("magnetic_flux_density", "value"),
        # Input('MAS_switch', 'value'),
        Input("isotope_id", "value"),
    ],
)
def update_plot(
    # submit_n_clicks
    # n1,
    spinning_frequency,
    number_of_points,
    frequency_bandwidth,
    reference_offset,
    magnetic_flux_density,
    # MAS_switch,
    isotope_id,
):
    """
    The method creates a new spectrum dictionary based on the inputs
    and re-computes the NMR lineshape.
    """

    if frequency_bandwidth in [None, 0, "", ".", "-"]:
        return [empty_plot(), "#ff2b2b"]
    if reference_offset in [None, "", ".", "-"]:
        return [empty_plot(), "#ff2b2b"]
    if spinning_frequency in [None, "", ".", "-"]:
        return [empty_plot(), "#ff2b2b"]

    # exit when the following conditions are True
    if number_of_points == 0 or isotope_id in ["", None]:
        return [empty_plot(), "#ff2b2b"]

    # calculating frequency_bandwidth
    try:
        frequency_bandwidth = float(frequency_bandwidth)
    except ValueError:
        return [empty_plot(), "#ff2b2b"]

    # calculating spin_frequency
    try:
        spin_frequency = float(spinning_frequency)
    except ValueError:
        return [empty_plot(), "#ff2b2b"]

    # calculating reference_offset
    try:
        reference_offset = float(reference_offset)
    except ValueError:
        return [empty_plot(), "#ff2b2b"]

    try:
        magnetic_flux_density = float(magnetic_flux_density) / 42.57747892
    except ValueError:
        return [empty_plot(), "#ff2b2b"]

    # if MAS_switch:
    rotor_angle_in_degree = 54.735

    spectrum = {
        "direct_dimension": {
            "isotope": isotope_id,
            "magnetic_flux_density": str(magnetic_flux_density * 100) + " T",
            "rotor_frequency": str(spin_frequency) + " kHz",
            "rotor_angle": str(rotor_angle_in_degree) + " deg",
            "number_of_points": 2 ** number_of_points,
            "spectral_width": str(frequency_bandwidth) + " kHz",
            "reference_offset": str(reference_offset) + " kHz",
        }
    }
    # print(spectrum["direct_dimension"])

    sim.spectrum = Spectrum.parse_json_with_units(spectrum)
    sim.x, sim.y = sim.run(one_d_spectrum, geodesic_polyhedron_frequency=75)
    sim.x = sim.x.to("kHz").value
    data = go.Scatter(x=sim.x, y=sim.y, mode="lines", opacity=1.0, name=isotope_id)

    x_label = str(isotope_id + f" frequency / kHz")

    return [
        {
            "data": [data],
            "layout": go.Layout(
                xaxis={
                    "type": "linear",
                    "title": x_label,
                    "ticks": "outside",
                    "showline": True,
                    "autorange": True,
                },
                yaxis={
                    "type": "linear",
                    "title": "arbitrary unit",
                    "ticks": "outside",
                    "showline": False,
                    "autorange": True,
                },
                autosize=False,
                transition={"duration": 100, "easing": "sin-in-out"},
                margin={"l": 50, "b": 40, "t": 5, "r": 5},
                legend={"x": 0, "y": 1},
                hovermode="closest",
            ),
        },
        "#00cc96",
    ]


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


@app.callback(
    [
        Output("filename_dataset", "children"),
        Output("data_time", "children"),
        Output("error_message", "children"),
        Output("isotope_id", "options"),
        Output("isotope_id", "value"),
    ],
    [Input("upload_data", "contents")],
    [State("upload_data", "filename"), State("upload_data", "last_modified")],
)
def update_isotopomers(content, filename, date):
    """Update the isotopomers when a new file is imported."""
    # FIRST = False
    # print(FIRST)
    children, success = parse_contents(content, filename, date)

    if success:
        # isotope_list = [
        #     dbc.DropdownMenuItem(isotope) for isotope in sim.unique_isotopes
        # ]
        isotope_list = [
            {"label": site_iso, "value": site_iso} for site_iso in sim.unique_isotopes
        ]
        isotope = isotope_list[0]["value"]
    else:
        sim.isotopomers = []
        isotope_list = []
        isotope = None
    return [children[0], children[1], children[2], isotope_list, isotope]


# @app.callback(
#     Output("spectrum_id_info", "children"), [Input("magnetic_flux_density", "value")]
# )
# def update_magnetic_flux_density(value):
#     """Update the value of magnetic flux density."""
#     return "@ {0} MHz".format("{0:.2f}".format(42.57747892 * float(value)))


# @app.callback(
#     [
#          Output("frequency_bandwidth", "valid"),
#          Output("frequency_bandwidth", "invalid")
#     ],
#     [Input("frequency_bandwidth", "value")],
# )
# def check_frequency_bandwidth(value):
#     if is_physical_quantity(value):
#         value = string_to_quantity(value)
#         if value > 0.0:
#             return [True, False]
#     return [False, True]


# @app.callback(
#     [Output("reference_offset", "valid"), Output("reference_offset", "invalid")],
#     [Input("reference_offset", "value")],
# )
# def check_reference_offset(value):
#     if value.isnumeric():
#         return [True, False]
#     return [False, True]


# @app.callback(
#     [Output("number_of_points", "valid"), Output("number_of_points", "invalid")],
#     [Input("number_of_points", "value")],
# )
# def check_number_of_points(value):
#     if value > 0:
#         return [True, False]
#     return [False, True]


@app.callback(
    Output("modal", "is_open"),
    [Input("info", "n_clicks"), Input("close", "n_clicks")],
    [State("modal", "is_open")],
)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open


# @app.callback(Output("spectrum_id", "children"), [Input("isotope_id", "value")])
# def update_spectrum_title(value):
#     """Update the title of the plot."""
#     if value is None:
#         return "Spectrum"
#     return "{} spectrum".format(value)


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
@app.callback(
    Output("navbar-collapse", "is_open"),
    [Input("navbar-toggler", "n_clicks")],
    [State("navbar-collapse", "is_open")],
)
def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


def empty_plot():
    data = go.Scatter(x=[-1, 1], y=[0, 0], text="", mode="lines", opacity=1.0)
    return {
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
            # legend={'x': 0, 'y': 1},
            hovermode="closest",
        ),
    }


def parse_contents(contents, filename, date):
    try:
        if "json" in filename:
            content_string = contents.split(",")[1]
            decoded = base64.b64decode(content_string)
            parse = json.loads(str(decoded, encoding="UTF-8"))["isotopomers"]
            # print(parse)
            sim.isotopomers = [Isotopomer.parse_json_with_units(item) for item in parse]
            return (
                [
                    filename,
                    datetime.datetime.fromtimestamp(date),
                    "Select a JSON serialized isotopomers file.",
                ],
                True,
            )

        else:
            return (
                ["", "", "A JSON file with valid list of isotopomers is required."],
                False,
            )

    except Exception:
        if FIRST:
            return (["", "", "Select a JSON serialized isotopomers file."], False)
        else:
            return ["", "", "There was an error reading the file."], False


if __name__ == "__main__":
    app.run_server(debug=True)
