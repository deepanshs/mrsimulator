# -*- coding: utf-8 -*-
import base64
import io
import json
import os

# import dash_html_components as html
import numpy as np
import flask
import plotly.graph_objs as go

from copy import deepcopy
from dash.dependencies import Input
from dash.dependencies import Output
from dash.dependencies import State

from mrsimulator import Isotopomer, Dimension
from mrsimulator.methods import one_d_spectrum

import dash_bootstrap_components as dbc
import dash_html_components as html

from mrsimulator.app.app import app, sim
from mrsimulator.app import navbar, sidebar
from mrsimulator.app.widgets import main_body

# from mrsimulator.app.post_simulation import line_broadening


__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]

app.layout = dbc.Container(
    [
        navbar.navbar_top,
        dbc.Row(
            [
                dbc.Col(sidebar.sidebar, xs=12, sm=12, md=12, lg=12, xl=3),
                dbc.Col(
                    [html.Div(main_body), html.Div(id="isotopomer_computed_log")],
                    xs=12,
                    sm=12,
                    md=12,
                    lg=12,
                    xl=9,
                ),
            ]
        ),
        navbar.navbar_bottom,
    ],
    fluid=True,
    style={"max-width": "1400px"},
    className="flex-display",
)
server = app.server


FIRST = True
count = 0
total_sites = 0


# update the link to the downloadable serialized file.
@app.callback(
    [Output("download_csdm", "href")],
    [Input("nmr_spectrum", "figure")],
    [
        State("filename_dataset", "children"),
        State("isotope_id-0", "value"),
        State("decompose", "active"),
        State("rotor_frequency-0", "value"),
        State("number_of_points-0", "value"),
        State("spectral_width-0", "value"),
        State("reference_offset-0", "value"),
        State("spectrometer_frequency-0", "value"),
    ],
)
def update_download_link(
    figure,
    filename_dataset,
    isotope_id,
    decompose_status,
    rotor_frequency,
    number_of_points,
    spectral_width,
    reference_offset,
    spectrometer_frequency,
):
    """Update the link to the downloadable serialized file."""
    name = os.path.splitext(str(filename_dataset))[0]

    global count
    count += 1

    head = "/dash/urlToDownload?name={0}&isotope={1}".format(name, isotope_id)
    tail = "&status={0}&index={1}_{2}_{3}_{4}_{5}_{6}".format(
        str(decompose_status),
        str(rotor_frequency),
        str(number_of_points),
        str(spectral_width),
        str(reference_offset),
        str(spectrometer_frequency),
        str(count),
    )
    return [head + tail]


# Serialize the computed spectrum and download the serialized file.
@app.server.route("/dash/urlToDownload")
def download():
    """Serialize the computed spectrum and download the serialized file."""
    # creating a dynamic csv or file here using `StringIO`
    filename = flask.request.args.get("name")
    isotope_id = flask.request.args.get("isotope")
    index = flask.request.args.get("index")
    is_decomposed = flask.request.args.get("status")

    str_io = io.StringIO()

    global count

    mem = io.BytesIO()
    file_name = "_".join([filename, isotope_id, is_decomposed, index, str(count)])

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
    [Output("nmr_spectrum", "figure")],
    [
        Input("rotor_frequency-0", "value"),
        Input("rotor_angle-0", "value"),
        Input("number_of_points-0", "value"),
        Input("spectral_width-0", "value"),
        Input("reference_offset-0", "value"),
        Input("spectrometer_frequency-0", "value"),
        Input("isotope_id-0", "value"),
        Input("decompose", "active"),
        # Input("broadening_points-0", "value"),
        Input("close_setting", "n_clicks"),
    ],
    [State("averaging_quality", "value")],
)
def update_data(
    rotor_frequency,
    rotor_angle,
    number_of_points,
    spectral_width,
    reference_offset,
    spectrometer_frequency,
    isotope_id,
    decompose,
    # broadening,
    close_setting_model_01,
    averaging_quality,
):
    """Evaluate the spectrum and update the plot."""
    clear_array = [clear_1D_plot()]

    if spectral_width in [None, 0, "", ".", "-"]:
        return clear_array
    if reference_offset in [None, "", ".", "-"]:
        return clear_array
    if rotor_frequency in [None, "", ".", "-"]:
        return clear_array

    # exit when the following conditions are True
    if number_of_points == 0 or isotope_id in ["", None]:
        return clear_array

    # calculating spectral_width
    try:
        spectral_width = float(spectral_width)
    except ValueError:
        return clear_array

    # calculating rotor_frequency
    try:
        rotor_frequency = float(eval(str(rotor_frequency)))
    except ValueError:
        return clear_array
    except SyntaxError:
        try:
            rotor_frequency = float(rotor_frequency)
        except ValueError:
            return clear_array

    # calculating reference_offset
    try:
        reference_offset = float(reference_offset)
    except ValueError:
        return clear_array

    try:
        magnetic_flux_density = float(spectrometer_frequency) / 42.57747892
    except ValueError:
        return clear_array

    # calculating rotor angle
    try:
        rotor_angle = float(eval(str(rotor_angle)))  # 54.735
    except ValueError:
        return clear_array
    except SyntaxError:
        return clear_array

    sim.nt = averaging_quality
    spectrum = {
        "isotope": isotope_id,
        "magnetic_flux_density": str(magnetic_flux_density * 100) + " T",
        "rotor_frequency": str(rotor_frequency) + " kHz",
        "rotor_angle": str(rotor_angle) + " deg",
        "number_of_points": 2 ** number_of_points,
        "spectral_width": str(spectral_width) + " kHz",
        "reference_offset": str(reference_offset) + " kHz",
        "nt": averaging_quality,
    }

    if not sim.spectrum_previous == spectrum or sim.new:
        sim.new = False
        sim.spectrum_previous = deepcopy(spectrum)

        # sim.spectrum = [Dimension()]
        sim.spectrum[0].isotope = isotope_id
        sim.spectrum[0].magnetic_flux_density = magnetic_flux_density * 100
        sim.spectrum[0].rotor_frequency = rotor_frequency * 1000.0
        sim.spectrum[0].rotor_angle = rotor_angle * np.pi / 180.0
        sim.spectrum[0].number_of_points = 2 ** number_of_points
        sim.spectrum[0].spectral_width = spectral_width * 1000.0
        sim.spectrum[0].reference_offset = reference_offset * 1000.0

        sim.x, sim.original = sim.run(
            one_d_spectrum,
            geodesic_polyhedron_frequency=sim.nt,
            individual_spectrum=True,
        )

        sim.x = sim.x.value
        sim.y = sim.original
    # remove the following line and add post simulation function as
    # sim.y = post_simulation(line_broadening, sigma=broadening)

    return plot_1D(isotope_id, decompose)


def post_simulation(function, **kwargs):
    return [
        function(datum, **kwargs)
        for datum in sim.original
        if not isinstance(datum, list)
    ]


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


@app.callback(
    Output("modal_setting", "is_open"),
    [Input("advance_setting", "n_clicks"), Input("close_setting", "n_clicks")],
    [State("modal_setting", "is_open")],
)
def toggle_modal_setting(n1, n2, is_open):
    """Model window for advance input."""
    if n1 or n2:
        return not is_open
    return is_open


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
    return [data_object]


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
            sim.spectrum = [Dimension(number_of_points=2, spectral_width=25)]
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
