# -*- coding: utf-8 -*-
import base64
import datetime
import io
import json
import os

import dash_core_components as dcc
import dash_html_components as html
import flask
import numpy as np
import plotly.graph_objs as go
from dash import Dash
from dash.dependencies import Input
from dash.dependencies import Output
from dash.dependencies import State

from mrsimulator import Simulator
from mrsimulator.methods import one_d_spectrum
from mrsimulator.widgets import direct_dimension_setup
from mrsimulator.widgets import display_isotopomers
from mrsimulator.widgets import plot_object_widget
from mrsimulator.widgets import spectrum_object_widget
from mrsimulator.widgets import top_bar


__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


external_stylesheets = [
    (
        "https://cdnjs.cloudflare.com/ajax/libs/materialize/1.0.0/css/"
        "materialize.min.css"
    )
]

external_scripts = [
    (
        "https://cdnjs.cloudflare.com/ajax/libs/materialize/1.0.0/js/"
        "materialize.min.js"
    )
]
colors = {"background": "#e2e2e2", "text": "#585858"}

__title__ = "mrsimulator"
__sub_title__ = "A web application framework for NMR lineshape simulation."


class mrsimulatorWebApp:
    __slots__ = "_filename"

    def __init__(self, app, simulator):
        self._filename = ""
        app.layout = html.Div(
            className="container",
            # style={'background-color': '#303030', 'color': 'white'},
            # style={"width": "90%"},
            children=[
                dcc.ConfirmDialog(
                    id="confirm", message="Cannot process the current request."
                ),
                html.H1(
                    children=__title__,
                    style={"textAlign": "center", "color": colors["text"]},
                ),
                html.Div(
                    children=__sub_title__,
                    style={"textAlign": "center", "color": colors["text"]},
                ),
                html.Hr(),
                *top_bar(),
                html.Hr(),
                html.Div(
                    className="row",
                    children=[
                        html.Div(
                            className="col s12 m12 l7", children=plot_object_widget()
                        ),
                        html.Div(
                            className="col s12 m12 l5",
                            children=spectrum_object_widget(direct_dimension_setup()),
                        ),
                    ],
                    # style={"width": "100%", "height": "100%"},
                ),
                html.Hr(),
                html.Div(id="isotopomer_computed_log"),
            ],
        )

        @app.callback(
            Output("download_link", "href"),
            [Input("nmr_spectrum", "figure")],
            [State("filename_dataset", "children"), State("spectrum_id", "children")],
        )
        def _update_link(value, filename_dataset, spectrum_id):
            """Update the csv download link when the plot is refreshed."""
            name_ = os.path.splitext(str(filename_dataset))[0]
            return "/dash/urlToDownload?value={0}+++{1}".format(name_, spectrum_id)

        @app.server.route("/dash/urlToDownload")
        def _download_csv():
            value = flask.request.args.get("value")

            # creating a dynamic csv or file here using `StringIO`
            str_io = io.StringIO()
            writer = np.asarray([simulator._freq.value, simulator._amp]).T
            file_, nuclei = value.split("   ")
            _header_ = (
                ("\n {0} from file {1}.json\nfrequency / kHz, amplitudes\n")
            ).format(nuclei, file_)

            # save file as csv
            np.savetxt(str_io, writer, fmt="%f", delimiter=",", header=_header_)

            mem = io.BytesIO()
            mem.write(str_io.getvalue().encode("utf-8"))
            mem.seek(0)
            str_io.close()
            name_ = "_".join([file_, nuclei])
            return flask.send_file(
                mem,
                mimetype="text/csv",
                attachment_filename=f"{name_}.csv",
                as_attachment=True,
            )

        @app.callback(
            Output("isotopomer_computed_log", "children"),
            [Input("nucleus_id", "value")],
        )
        def _update_isotopomers_table_log(value):
            return display_isotopomers(value, simulator.isotopomers)

        # @app.callback(
        #     Output('tabs_content', 'children'),
        #     [Input('tabs', 'value')]
        # )
        # def update_tab(value):
        #     if value == 'direct_dimension':
        #         return direct_dimension_setup()

        @app.callback(
            [Output("nmr_spectrum", "figure")],
            [
                # Input('confirm', 'submit_n_clicks'),
                Input("spinning_frequency_in_kHz_coarse", "value"),
                Input("spinning_frequency_in_kHz_fine", "value"),
                Input("number_of_points", "value"),
                Input("frequency_bandwidth_coarse", "value"),
                Input("frequency_bandwidth_fine", "value"),
                Input("reference_offset_coarse", "value"),
                Input("reference_offset_fine", "value"),
                Input("magnetic_flux_density", "value"),
                # Input('MAS_switch', 'value'),
                Input("nucleus_id", "value"),
            ],
        )
        def _update_plot(
            # submit_n_clicks,
            spinning_frequency_in_kHz_coarse,
            spinning_frequency_in_kHz_fine,
            number_of_points,
            frequency_bandwidth_coarse,
            frequency_bandwidth_fine,
            reference_offset_coarse,
            reference_offset_fine,
            magnetic_flux_density,
            # MAS_switch,
            nucleus_id,
        ):
            """
            The method creates a new spectrum dictionary based on the inputs
            and re-computes the NMR lineshape.
            """

            # calculating frequency_bandwidth
            frequency_bandwidth = frequency_bandwidth_coarse + frequency_bandwidth_fine

            # exit when the following conditions are True
            if (
                number_of_points == 0
                or frequency_bandwidth == 0
                or nucleus_id in ["", None]
            ):
                return empty_plot()

            # calculating spin_frequency
            spin_frequency = (
                spinning_frequency_in_kHz_coarse + spinning_frequency_in_kHz_fine
            )

            reference_offset = reference_offset_coarse + reference_offset_fine

            # if MAS_switch:
            rotor_angle_in_degree = 54.735
            # else:
            #     rotor_angle_in_degree = 0

            simulator.spectrum = {
                "direct_dimension": {
                    "nucleus": nucleus_id,
                    "magnetic_flux_density": str(magnetic_flux_density) + " T",
                    "rotor_frequency": str(spin_frequency) + " kHz",
                    "rotor_angle": str(rotor_angle_in_degree) + " deg",
                    "number_of_points": 2 ** number_of_points,
                    "spectral_width": str(frequency_bandwidth) + " kHz",
                    "reference_offset": str(reference_offset) + " kHz",
                }
            }

            freq, amp = simulator.run(one_d_spectrum)
            freq = freq.to("kHz")
            data_spinning = go.Scatter(
                x=freq, y=amp / amp.max(), mode="lines", opacity=1.0, name=nucleus_id
            )

            x_label = str(nucleus_id + f" frequency / {freq.unit}")

            return [
                {
                    "data": [data_spinning],
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
                            "showline": True,
                            "autorange": True,
                        },
                        autosize=True,
                        # 'easing': 'quad-in-out'},
                        transition={"duration": 500},
                        margin={"l": 50, "b": 40, "t": 5, "r": 5},
                        # legend={'x': 0, 'y': 1},
                        hovermode="closest",
                    ),
                }
            ]

        @app.callback(
            [
                Output("filename_dataset", "children"),
                Output("data_time", "children"),
                Output("error_message", "children"),
                Output("nucleus_widget_id", "children"),
            ],
            [Input("upload_data", "contents")],
            [State("upload_data", "filename"), State("upload_data", "last_modified")],
        )
        def _update_isotopomers(content, filename, date):
            """Update the isotopomers when a new file is imported."""
            # FIRST = False
            children, success = parse_contents(simulator, content, filename, date)

            if success:
                nuclei = [
                    {"label": site_iso, "value": site_iso}
                    for site_iso in simulator.isotope_list
                ]

                if len(simulator.isotope_list) >= 1:
                    value = simulator.isotope_list[0]
                else:
                    value = ""
                nucleus_id = [
                    dcc.Dropdown(
                        id="nucleus_id",
                        options=nuclei,
                        value=value,
                        style={
                            "display": "block",
                            "margin-left": "auto",
                            "margin-right": "auto",
                            "width": "auto",
                        },
                    )
                ]
            else:
                simulator.isotopomers = []
                nucleus_id = [
                    dcc.Dropdown(
                        id="nucleus_id",
                        style={
                            "display": "block",
                            "margin-left": "auto",
                            "margin-right": "auto",
                            "width": "auto",
                        },
                    )
                ]
            return [children[0], children[1], children[2], nucleus_id]

        @app.callback(
            Output("Magnetic_flux_density_output_container", "children"),
            [Input("magnetic_flux_density", "value")],
        )
        def _update_magnetic_flux_density(value):
            """Update the value of magnetic flux density."""
            return "Magnetic flux density   {0} T @ {1} MHz".format(
                value, "{0:.2f}".format(42.57747892 * value)
            )

        @app.callback(Output("spectrum_id", "children"), [Input("nucleus_id", "value")])
        def _update_spectrum_title(value):
            """Update the title of the plot."""
            if value is None:
                return "Spectrum"
            return "{} spectrum".format(value)

        @app.callback(
            Output("number_of_points_output_container", "children"),
            [Input("number_of_points", "value")],
        )
        def _update_number_of_points(value):
            """Update the number of points."""
            return "Number of points        {}".format(2 ** value)

        @app.callback(
            Output("spinning_frequency_output_container", "children"),
            [
                Input("spinning_frequency_in_kHz_fine", "value"),
                Input("spinning_frequency_in_kHz_coarse", "value"),
            ],
        )
        def _update_rotor_frequency(value1, value2):
            """Update the rotor spin frequency."""
            return "Magic angle spinning frequency {} kHz".format(value1 + value2)

        @app.callback(
            Output("reference_offset_output_container", "children"),
            [
                Input("reference_offset_fine", "value"),
                Input("reference_offset_coarse", "value"),
            ],
        )
        def _update_reference_offset(value1, value2):
            """Update the reference offset."""
            return "Reference offset {} kHz".format(value1 + value2)

        @app.callback(
            Output("frequency_bandwidth_output_container", "children"),
            [
                Input("frequency_bandwidth_fine", "value"),
                Input("frequency_bandwidth_coarse", "value"),
            ],
        )
        def _update_frequency_bandwidth(value1, value2):
            """Update the spectral width."""
            return "Spectral width {} kHz".format(value1 + value2)


#  ['linear', 'quad', 'cubic', 'sin', 'exp', 'circle',
#             'elastic', 'back', 'bounce', 'linear-in', 'quad-in',
#             'cubic-in', 'sin-in', 'exp-in', 'circle-in', 'elastic-in',
#             'back-in', 'bounce-in', 'linear-out', 'quad-out',
#             'cubic-out', 'sin-out', 'exp-out', 'circle-out',
#             'elastic-out', 'back-out', 'bounce-out', 'linear-in-out',
#             'quad-in-out', 'cubic-in-out', 'sin-in-out', 'exp-in-out',
#             'circle-in-out', 'elastic-in-out', 'back-in-out',
#             'bounce-in-out']


def empty_plot():
    data = go.Scatter(x=[-1, 1], y=[0, 0], text="", mode="lines", opacity=1.0)
    return [
        {
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
    ]


def parse_contents(simulator, contents, filename, date):
    try:
        if "json" in filename:
            content_string = contents.split(",")[1]
            decoded = base64.b64decode(content_string)
            parse = json.loads(str(decoded, encoding="UTF-8"))["isotopomers"]
            # print(parse)
            simulator.isotopomers = parse
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
    app = Dash(__name__)
    for css in external_stylesheets:
        app.css.append_css({"external_url": css})
    for js in external_scripts:
        app.scripts.append_script({"external_url": js})

    FIRST = True
    simulator = Simulator()
    mrsimulatorWebApp(app, simulator)
    app.run_server(debug=True)
