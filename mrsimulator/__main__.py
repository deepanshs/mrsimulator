# -*- coding: utf-8 -*-

from mrsimulator import Simulator
from mrsimulator.methods import one_d_spectrum
from dash import Dash
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
import plotly.graph_objs as go
from dash.dependencies import Input, Output, State
import datetime, json
import base64

from mrsimulator.widgets import (
    spectrum_object_widget,
    plot_object_widget,
    direct_dimension_setup,
)


__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.90@osu.edu"


external_stylesheets = [
    "https://cdnjs.cloudflare.com/ajax/libs/materialize/1.0.0/css/materialize.min.css"
]


external_scripts = [
    'https://cdnjs.cloudflare.com/ajax/libs/materialize/1.0.0/js/materialize.min.js'
]


colors = {
    'background': '#e2e2e2',
    'text': '#585858'
}



# isotopomer = get_system()
def setup(app, simulator):
    app.layout = html.Div(className='container',
        # style={'background-color': '#303030', 'color': 'white'},
        # style={'backgroundColor': colors['background']},
        children=[
            # dcc.ConfirmDialog(
            #     id='confirm',
            #     message='Cannot process the current requeste.',
            # ),
            html.H1(
                children='mrsimulator',
                style={
                    'textAlign': 'center',
                    'color': colors['text']
                }
            ),
            html.Div(children='A web application framework for NMR lineshape simulation.',
                style={
                    'textAlign': 'center',
                    'color': colors['text']
                }
            ),
            html.Hr(),
            html.Div(className='row', children=[
                html.Div(className='col s12 m7 l7', id='output-data-upload'),
                html.Div(className='col s12 m5 l5', children=[
                    dcc.Upload(
                        id='upload_data',
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select File')
                        ]),
                        style={
                            'width': '100%',
                            'height': '50px',
                            'lineHeight': '50px',
                            'borderWidth': '1px',
                            'borderStyle': 'dashed',
                            'borderRadius': '5px',
                            'textAlign': 'center',
                            'margin': '1px'
                        },
                        # Allow multiple files to be uploaded
                        multiple=False
                    ),
                    html.Div(id='error_message', style={'textAlign': 'center'})
                ])
            ]),
            html.Hr(),
            html.Div(
                className='row',
                children=[
                    html.Div(
                        className='col s12 m7 l7',
                        children=plot_object_widget()
                    ),
                    html.Div(
                        className='col s12 m5 l5',
                        children=spectrum_object_widget(direct_dimension_setup())#nuclei, simulator.isotope_list)
                    )
            ]),
        ]
    )

    # @app.callback(
    #     Output('tabs_content', 'children'),
    #     [Input('tabs', 'value')]
    # )
    # def update_tab(value):
    #     if value == 'direct_dimension':
    #         return direct_dimension_setup()

    @app.callback(
        [
            Output('nmr_spectrum', 'figure'),
        ],
        [
            # Input('confirm', 'submit_n_clicks'),
            Input(component_id='spinning_frequency_in_kHz_coarse', component_property='value'),
            Input(component_id='spinning_frequency_in_kHz_fine', component_property='value'),
            Input(component_id='number_of_points', component_property='value'),
            Input(component_id='frequency_bandwidth_coarse', component_property='value'),
            Input(component_id='frequency_bandwidth_fine', component_property='value'),
            Input(component_id='reference_offset_coarse', component_property='value'),
            Input(component_id='reference_offset_fine', component_property='value'),
            Input(component_id='magnetic_flux_density', component_property='value'),
            # Input(component_id='MAS_switch', component_property='value'),
            Input(component_id='nucleus_id', component_property='value'),
        ])
    def update_plot(
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
            nucleus_id):

        frequency_bandwidth = frequency_bandwidth_coarse + frequency_bandwidth_fine
        if number_of_points == 0 or frequency_bandwidth == 0 or nucleus_id in ['', None]:
            return empty_plot()

        spin_frequency = spinning_frequency_in_kHz_coarse+spinning_frequency_in_kHz_fine
        reference_offset = reference_offset_coarse + reference_offset_fine
        
        # if MAS_switch:
        rotor_angle_in_degree = 54.735
        # else:
        #     rotor_angle_in_degree = 0

        simulator.spectrum = {
            "direct_dimension": {
                "nucleus": nucleus_id,
                "magnetic_flux_density": str(magnetic_flux_density) + " T",
                "rotor_frequency": str(spin_frequency)+' kHz',
                "rotor_angle": str(rotor_angle_in_degree)+' deg',
                "number_of_points": 2**number_of_points,
                "spectral_width": str(frequency_bandwidth)+' kHz',
                "reference_offset": str(reference_offset)+' kHz',
            }
        }

        freq, amp = simulator.run(one_d_spectrum)

        data_spinning = go.Scatter(x=freq/1000, y=amp/amp.max(),
                                mode='lines', opacity=1.0, name=nucleus_id)
        x_label = str(nucleus_id + ' frequency / kHz')
        # print(x_label)
        return [{
            'data': [data_spinning],
            'layout': go.Layout(
                xaxis={
                    'type': 'linear',
                    'title': x_label, #'frequency / kHz',
                    'ticks': 'outside',
                    'showline': True,
                    'autorange': True
                },
                yaxis={
                    'type': 'linear',
                    'title': 'arbitrary unit',
                    'ticks': 'outside',
                    'showline': True,
                    'autorange': True,
                },
                autosize=True,
                transition={'duration': 600}, #'easing': 'quad-in-out'},
                margin={'l': 50, 'b': 40, 't': 5, 'r': 5},
                # legend={'x': 0, 'y': 1},
                hovermode='closest'
                )
            }]

    @app.callback(
        [
            Output('output-data-upload', 'children'),
            Output('error_message', 'children'),
            Output('nucleus_widget_id', 'children'),
        ],
        [Input('upload_data', 'contents')],
        [State('upload_data', 'filename'),
        State('upload_data', 'last_modified')])
    def update_isotopomers(list_of_contents, list_of_names, list_of_dates):
        # if list_of_contents is not None:
        FIRST = False
        children, success = parse_contents(simulator, list_of_contents, 
                            list_of_names,
                            list_of_dates)
        

        if success:
            nuclei =[
                {'label': site_iso, 'value': site_iso}
                    for site_iso in simulator.isotope_list
            ]
            
            if len(simulator.isotope_list) >= 1:
                value = simulator.isotope_list[0]
            else:
                value = ''
            nucleus_id = [
                    dcc.Dropdown(
                        id='nucleus_id',
                        options=nuclei,
                        value=value,
                        style={"display": "block",
                            "margin-left": "auto",
                            "margin-right": "auto",
                            "width": "auto"},
                    ),
                ]
        else:
            simulator.isotopomers = []
            nucleus_id = [
                    dcc.Dropdown(
                        id='nucleus_id',
                        style={"display": "block",
                            "margin-left": "auto",
                            "margin-right": "auto",
                            "width": "auto"},
                    ),
                ]
        return [children[0]], [children[1]], nucleus_id

    @app.callback(
        Output('Magnetic_flux_density_output_container', 'children'),
        [Input('magnetic_flux_density', 'value')]
    )
    def update_magnetic_flux_density(value):
        return 'Magnetic flux density   {0} T @ {1} MHz'.format(
                            value, '{0:.2f}'.format(42.57747892*value)
                )

    @app.callback(
        Output('spectrum_id', 'children'),
        [Input('nucleus_id', 'value')]
    )
    def update_spectrum_title(value):
        if value is None:
            return 'Spectrum'
        return '{} spectrum'.format(value)

    @app.callback(
        Output('number_of_points_output_container', 'children'),
        [Input('number_of_points', 'value')]
    )
    def update_number_of_points(value):
        return 'Number of points        {}'.format(2**value)

    @app.callback(
        Output('spinning_frequency_output_container', 'children'),
        [Input('spinning_frequency_in_kHz_fine', 'value'),
         Input('spinning_frequency_in_kHz_coarse', 'value')]
    )
    def update_rotor_frequency(value1, value2):
        return 'Magic angle spinning    {} kHz'.format(value1 + value2)

    @app.callback(
        Output('reference_offset_output_container', 'children'),
        [Input('reference_offset_fine', 'value'),
         Input('reference_offset_coarse', 'value')]
    )
    def update_reference_offset(value1, value2):
        return 'Reference offset        {} kHz'.format(value1 + value2)

    @app.callback(
        Output('frequency_bandwidth_output_container', 'children'),
        [Input('frequency_bandwidth_fine', 'value'),
         Input('frequency_bandwidth_coarse', 'value')]
    )
    def update_frequency_bandwidth(value1, value2):
        return 'Spectral width          {} kHz'.format(value1 + value2)


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
    data = go.Scatter(x=[-1,1], y=[0,0], text='', mode='lines', opacity=1.0)
    # return [{'data': [data]}]
    return [{
        'data': [data],
        'layout': go.Layout(
            xaxis={
                'type': 'linear',
                'title': 'frequency / kHz',
                'ticks': 'outside',
                'showline': True,
                'autorange': True
            },
            yaxis={
                'type': 'linear',
                'title': 'arbitrary unit',
                'ticks': 'outside',
                'showline': True,
                'autorange': True,
            },
            autosize=True,
            transition={'duration': 500},
            margin={'l': 50, 'b': 40, 't': 5, 'r': 5},
            # legend={'x': 0, 'y': 1},
            hovermode='closest'
            )
        }]


def parse_contents(simulator, contents, filename, date):
    try:
        if 'json' in filename:
            content_type, content_string = contents.split(',')
            decoded = base64.b64decode(content_string)
            parse=json.loads(str(decoded, encoding="UTF-8"))['isotopomers']
            # print(parse)
            simulator.isotopomers = parse
            return [html.Div([
                html.H5(filename),
                html.H6(datetime.datetime.fromtimestamp(date))
            ]), html.Label([
                'Select a JSON serialized isotopomers file.'
            ])], True

        else:
            return ['', html.Label([
                'A JSON file is required.'
            ])], False
    except: #Exception as e:
        # simulator.isotopomers = []
        # print(e)
        if FIRST:
            return ['', html.Label([
                'Select a JSON serialized isotopomers file.'
            ])], False
        else:
            return ['', html.Label([
                'There was an error reading the file.'
            ])], False




if __name__ == '__main__':
    app = Dash(__name__)
    for css in external_stylesheets:
        app.css.append_css({"external_url": css})
    for js in external_scripts:
        app.scripts.append_script({'external_url': js})
    
    FIRST = True
    simulator = Simulator()
    setup(app, simulator)
    app.run_server(debug=True)
