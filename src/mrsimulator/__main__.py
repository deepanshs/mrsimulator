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

from mrsimulator import Simulator, Isotopomer, SpectroscopicDimension
from mrsimulator.methods import one_d_spectrum
from mrsimulator.web_ui.widgets import get_isotopomers, input_file, main_body
from mrsimulator.unit import string_to_quantity, is_physical_quantity

import csdmpy as cp

from mrsimulator.web_ui import navbar, titlebar


__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


external_scripts = [
    "https://cdnjs.cloudflare.com/ajax/libs/materialize/1.0.0/js/materialize.min.js"
]


app = Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    external_scripts=external_scripts,
)


FIRST = False
sim = Simulator()
sim.x = np.asarray([-1.0, 1.0])
sim.y = np.asarray([0.0, 0.0])
sim.nt = 64
count = 0
total_sites = 0

app.layout = dbc.Container(
    [
        navbar.navbar,
        html.Br(),
        html.Div(input_file),
        html.Br(),
        html.Div(main_body, className="accordion"),
        html.Br(),
        # html.Div(id="isotopomer_computed_log"),
        dbc.Jumbotron(),
    ],
    style={
        "color": "#585858",
        "background-image": "linear-gradient(#fdfefe, #f4f6f7, #ffffff)",
    },
)


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
        State("decompose", "active"),
        State("spinning_frequency", "value"),
        State("number_of_points", "value"),
        State("frequency_bandwidth", "value"),
        State("reference_offset", "value"),
        State("magnetic_flux_density", "value"),
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
    """Update the csv download link when the plot is refreshed."""
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


@app.server.route("/dash/urlToDownload")
def download():
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
        new = cp.new()
        dimension = {
            "type": "linear",
            "count": sim.spectrum[0].number_of_points,
            "increment": "{0} Hz".format(
                sim.spectrum[0].spectral_width / sim.spectrum[0].number_of_points
            ),
            "coordinates_offset": f"{sim.spectrum[0].reference_offset} Hz",
            "origin_offset": f"{sim.spectrum[0].larmor_frequency} MHz",
            "complex_fft": True,
        }
        new.add_dimension(dimension)

        if is_decomposed == "True":
            for index, datum in enumerate(sim.y):
                if not isinstance(datum, list):

                    dependent_variable = {
                        "type": "internal",
                        "quantity_type": "scalar",
                        "numeric_type": "float64",
                        "components": [datum],
                    }

                    name = sim.isotopomers[index].name
                    if name != "":
                        dependent_variable.update({"name": sim.isotopomers[index].name})

                    description = sim.isotopomers[index].description
                    if description != "":
                        dependent_variable.update(
                            {"description": sim.isotopomers[index].description}
                        )
                    new.add_dependent_variable(dependent_variable)
                    new.dependent_variables[-1].encoding = "base64"

        else:
            sum_spectrums = np.zeros(sim.spectrum[0].number_of_points)
            for index, datum in enumerate(sim.y):
                if not isinstance(datum, list):
                    sum_spectrums += datum

            dependent_variable = {
                "type": "internal",
                "quantity_type": "scalar",
                "numeric_type": "float64",
                "components": [sum_spectrums],
            }
            new.add_dependent_variable(dependent_variable)
            new.dependent_variables[-1].encoding = "base64"

        new.save(output_device=str_io)

        mem.write(str_io.getvalue().encode("utf-8"))
        mem.seek(0)
        str_io.close()
        return flask.send_file(
            mem, attachment_filename=f"{file_name}.csdf", as_attachment=True
        )


@app.callback(
    [Output("nmr_spectrum", "figure"), Output("indicator_status", "color")],
    [
        Input("spinning_frequency", "value"),
        Input("number_of_points", "value"),
        Input("frequency_bandwidth", "value"),
        Input("reference_offset", "value"),
        Input("magnetic_flux_density", "value"),
        Input("isotope_id", "value"),
        Input("decompose", "active"),
    ],
)
def update_plot(
    spinning_frequency,
    number_of_points,
    frequency_bandwidth,
    reference_offset,
    magnetic_flux_density,
    isotope_id,
    show_individual,
):
    """
    The method creates a new spectrum dictionary based on the inputs
    and re-computes the NMR lineshape.
    """
    color = "danger"  # "#ff2b2b"
    if frequency_bandwidth in [None, 0, "", ".", "-"]:
        return [empty_plot(), color]
    if reference_offset in [None, "", ".", "-"]:
        return [empty_plot(), color]
    if spinning_frequency in [None, "", ".", "-"]:
        return [empty_plot(), color]

    # exit when the following conditions are True
    if number_of_points == 0 or isotope_id in ["", None]:
        return [empty_plot(), color]

    # calculating frequency_bandwidth
    try:
        frequency_bandwidth = float(frequency_bandwidth)
    except ValueError:
        return [empty_plot(), color]

    # calculating spin_frequency
    try:
        spin_frequency = float(spinning_frequency)
    except ValueError:
        return [empty_plot(), color]

    # calculating reference_offset
    try:
        reference_offset = float(reference_offset)
    except ValueError:
        return [empty_plot(), color]

    try:
        magnetic_flux_density = float(magnetic_flux_density) / 42.57747892
    except ValueError:
        return [empty_plot(), color]

    # if MAS_switch:
    rotor_angle_in_degree = 54.735

    spectrum = {
        "isotope": isotope_id,
        "magnetic_flux_density": str(magnetic_flux_density * 100) + " T",
        "rotor_frequency": str(spin_frequency) + " kHz",
        "rotor_angle": str(rotor_angle_in_degree) + " deg",
        "number_of_points": 2 ** number_of_points,
        "spectral_width": str(frequency_bandwidth) + " kHz",
        "reference_offset": str(reference_offset) + " kHz",
    }

    sim.spectrum = [SpectroscopicDimension.parse_dict_with_units(spectrum)]
    sim.x, sim.y = sim.run(
        one_d_spectrum, geodesic_polyhedron_frequency=sim.nt, individual_spectrum=True
    )

    sim.x = sim.x.value

    # print(sim.y)

    data = []
    if show_individual:
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
                        # name=isotope_id,
                        name=name,
                        # hovername=
                        # customdata=str(sim.isotopomers[i]),
                    )
                )
    else:
        sum_spectrums = np.zeros(sim.spectrum[0].number_of_points)
        for datum in sim.y:
            if not isinstance(datum, list):
                sum_spectrums += datum
        data.append(
            go.Scatter(
                x=sim.x,
                y=sum_spectrums,
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
            autosize=False,
            transition={"duration": 150, "easing": "sin-out"},
            margin={"l": 50, "b": 40, "t": 5, "r": 5},
            legend={"x": 0, "y": 1},
            # legend_orientation="h",
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
    children, success = parse_contents(content, filename, date)

    if success:
        isotope_list = [
            {"label": site_iso, "value": site_iso} for site_iso in sim.get_isotopes()
        ]
        isotope = isotope_list[0]["value"]
    else:
        sim.isotopomers = []
        isotope_list = []
        isotope = None
    return [children[0], children[1], children[2], isotope_list, isotope]


# @app.callback(
#     Output("modal", "is_open"),
#     [Input("info", "n_clicks"), Input("close", "n_clicks")],
#     [State("modal", "is_open")],
# )
# def toggle_modal(n1, n2, is_open):
#     if n1 or n2:
#         return not is_open
#     return is_open


@app.callback(
    Output("modal_setting", "is_open"),
    [Input("advance-id", "n_clicks"), Input("close_setting", "n_clicks")],
    [State("modal_setting", "is_open")],
)
def toggle_modal_setting(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open


@app.callback(Output("spectrum_id", "children"), [Input("isotope_id", "value")])
def update_spectrum_title(value):
    """Update the title of the plot."""
    if value is None:
        return "Spectrum"
    return "{} spectrum".format(value)


@app.callback(
    Output("number_of_averaging_points", "children"),
    [Input("averaging_quality", "value")],
)
def update_number_of_orientations(value):
    """Update the number of orientation for powder averaging."""
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
@app.callback(
    Output("navbar-collapse", "is_open"),
    [Input("navbar-toggler", "n_clicks")],
    [State("navbar-collapse", "is_open")],
)
def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


@app.callback(
    [
        Output(f"collapse-dimension", "is_open"),
        # *[Output(f"collapse-card-{i}", "n_clicks") for i in range(total_sites)],
    ],
    [
        Input(f"dimension-toggle", "n_clicks"),
        # *[Input(f"card-toggle-{i}", "n_clicks") for i in range(total_sites)],
    ],
    [
        State(f"collapse-dimension", "is_open"),
        # *[State(f"collapse-card-{i}", "n_clicks") for i in range(total_sites)],
    ],
)
def toggle_accordion(n1, is_dimension_open):
    ctx = dash.callback_context
    # global total_sites
    # print("total_sites", total_sites)
    # list_output = [False for _ in range(total_sites + 1)]

    if not ctx.triggered:
        return [""]
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if button_id == "dimension-toggle" and n1:
        # list_output[0] = not is_dimension_open
        return [not is_dimension_open]
    # for i in range(total_sites):
    #     if button_id == "card-toggle-{i}" and n_list[i]:
    #         list_output[i] = not is_open[i]
    #     return False, not is_open2, False
    # elif button_id == "group-3-toggle" and n3:
    #     return False, False, not is_open3
    return [False]


@app.callback(
    [Output("decompose", "active")],
    [Input("decompose", "n_clicks")],
    [State("decompose", "active")],
)
def update_spectral_decomposition(n1, status):
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


def empty_plot():
    data = go.Scatter(x=[-1, 1], y=[0, 0], text="", mode="lines", opacity=1.0)
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


def parse_contents(contents, filename, date):
    try:
        if "json" in filename:
            content_string = contents.split(",")[1]
            decoded = base64.b64decode(content_string)
            parse = json.loads(str(decoded, encoding="UTF-8"))["isotopomers"]
            sim.isotopomers = [Isotopomer.parse_dict_with_units(item) for item in parse]
            sim.spectrum = []
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
            return ["", "", "Error reading the file."], False


if __name__ == "__main__":
    app.run_server(debug=True)
