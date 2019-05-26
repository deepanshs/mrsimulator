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


def child_1():
    """
    Return the layout for nucleus, number of points,
    spectral width and reference offset.
    """
    return [
        html.Div(
            className='card-panel hoverable',
            children=[
                html.H5(
                    children='Spectrum parameters',
                    style={
                        'textAlign': 'left',
                        'color': colors['text']
                    }
                ),
                html.H6(
                    children='Direct dimensions',
                    style={
                        'textAlign': 'left',
                        'color': colors['text']
                    }
                ),
                html.Div(className='row', children=[
                    html.Div(className='col', children=[html.Label('Nucleus')]),
                    # html.Div(className='col', children=[daq.Indicator(
                    #     id='simulation_success_indicator',
                    #     size=10,
                    #     value=True,
                    #     color="#00cc96"
                    # )]),
                ]),
                html.Div(id='nucleus_widget_id', children=[
                    dcc.Dropdown(
                        id='nucleus_id',
                        style={"display": "block",
                            "margin-left": "auto",
                            "margin-right": "auto",
                            "width": "auto"},
                    ),
                ]),
                html.Div([
                    html.Div([
                        html.Label('Number of points'),
                    ]),
                    html.Div([
                        dcc.Input(id='number_of_points',
                                value=1024,
                                type='number',
                                min=0.0,
                                #   style={'width': '50%', 'display': 'inline-block'},
                        ),
                    ]),
                ]),
                html.Div([
                    html.Div([
                        html.Label('spectral width in Hz'),
                    ]),
                    html.Div([
                        dcc.Input(id='frequency_bandwidth',
                                value=100000.0,
                                type='number',
                                min=0.0,
                                #   style={'width': '50%', 'display': 'inline-block'},
                        ),
                    ]),
                ]),
                html.Div([
                    html.Div([
                        html.Label('reference offset in Hz'),
                    ]),
                    html.Div([
                        dcc.Input(id='reference_offset',
                                value=0,
                                type='number',
                                min=0.0,
                                #   style={'width': '50%', 'display': 'inline-block'},
                        ),
                    ]),
                ]),
                html.Div([
                    html.Div([
                        html.Label('Magic angle spinning'),
                    ]),
                    html.Div([
                        dcc.Slider(
                            id='spinning_frequency_in_kHz',
                            min=0.0,
                            max=50,
                            step=0.050,
                            value=0,
                        ),
                        html.Div(id='spinning_frequency_output_container'),
                    ]),
                ]),
            ]
        )
    ]


def plot_layout():
    return [
        html.Div(
            className='card-panel hoverable',
            children=[
                html.H5(
                    children='Spectrum',
                    style={
                        'textAlign': 'left',
                        'color': colors['text']
                    }
                ),
                dcc.Graph(
                    id='nmr_spectrum',
                    # animate=True,
                    # style={"display": "inline-block",
                    #     "margin-left": 'auto',
                    #     "margin-right": 'auto',
                    #     "width": 'auto'},
                    figure={ 'data': []}
                )
            ]
        )
    ]

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
            html.Div(children='A web application framework for NMR spectrum simulator.',
                style={
                    'textAlign': 'center',
                    'color': colors['text']
                }
            ),
            html.Hr(),
            html.Div(className='row', children=[
                html.Div(className='col s8 m8 l8', id='output-data-upload'),
                html.Div(className='col s4 m4 l4', children=[
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
                    html.Div(id='error_message')
                ])
            ]),
            html.Hr(),
            html.Div(
                className='row',
                children=[
                    html.Div(
                        className='col s12 m8 l8',
                        children=plot_layout()
                    ),
                    html.Div(
                        className='col s12 m4 l4',
                        children=child_1()#nuclei, simulator.isotope_list)
                    )
            ]),
        ]
    )

    @app.callback(
        [
            Output('nmr_spectrum', 'figure'),
        ],
        [
            # Input('confirm', 'submit_n_clicks'),
            Input(component_id='spinning_frequency_in_kHz', component_property='value'),
            Input(component_id='number_of_points', component_property='value'),
            Input(component_id='frequency_bandwidth', component_property='value'),
            Input(component_id='reference_offset', component_property='value'),
            # Input(component_id='MAS_switch', component_property='value'),
            Input(component_id='nucleus_id', component_property='value'),
        ])
    def update_plot(
            # submit_n_clicks,
            spinning_frequency_in_kHz,
            number_of_points,
            frequency_bandwidth,
            reference_offset,
            # MAS_switch,
            nucleus_id):

        if number_of_points == 0 or frequency_bandwidth == 0 or nucleus_id in ['', None]:
            return empty_plot()

        # if MAS_switch:
        rotor_angle_in_degree = 54.735
        # else:
        #     rotor_angle_in_degree = 0

        simulator.spectrum = {
            "direct_dimension": {
                "nucleus": nucleus_id,
                "magnetic_flux_density": "9.4 T",
                "rotor_frequency": str(spinning_frequency_in_kHz)+' kHz',
                "rotor_angle": str(rotor_angle_in_degree)+' deg',
                "number_of_points": number_of_points,
                "spectral_width": str(frequency_bandwidth)+' Hz',
                "reference_offset": str(reference_offset)+' Hz',
            }
        }

        freq, amp = simulator.run(one_d_spectrum)

        data_spinning = go.Scatter(x=freq/1000, y=amp/amp.max(),
                                mode='lines', opacity=1.0, name=nucleus_id)
        # x_label = str(nucleus_id + ' frequency / kHz')
        # print(x_label)
        return [{
            'data': [data_spinning],
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
                transition={'duration': 500, 'easing': 'quad-in-out'},
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
        Output('spinning_frequency_output_container', 'children'),
        [Input('spinning_frequency_in_kHz', 'value')]
    )
    def update_rotor_frequency(value):
        return html.Label('{} kHz'.format(value))


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
            transition={'duration': 500, 'easing': 'quad-in-out'},
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
    app.run_server()#debug=True)
