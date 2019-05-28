

import dash_html_components as html
import dash_core_components as dcc
import dash_daq as daq


__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.90@osu.edu"


colors = {
    'background': '#e2e2e2',
    'text': '#585858'
}


def plot_object_widget():
    return [
        html.Div(
            className='card-panel hoverable',
            children=[
                html.H4(
                    id='spectrum_id',
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
                ),
            ]
        )
    ]


def spectrum_object_widget(object_=[]):
    """
    Return the layout for nucleus, number of points,
    spectral width and reference offset.
    """        
    return [
        html.Div(
            className='card-panel hoverable',
            children=[
                html.H4(
                    children='Spectrum parameters',
                    style={
                        'textAlign': 'left',
                        'color': colors['text']
                    }
                ),

                # dcc.Tabs(id="tabs", value='direct_dimension', children=[
                #     dcc.Tab(label='Direct dimension', value='direct_dimension'),
                #     # dcc.Tab(label='Tab Two', value='tab-2-example'),
                # ]),
                html.Div(id='tabs_content', children=object_)
            ]
        )
    ]   


def direct_dimension_setup():
    return [
        html.H5(
            children='Direct dimension',
            style={
                'textAlign': 'left',
                'color': colors['text']
            }
        ),

        # environment
        html.Div([html.H6('Environment parameters')],
            style={
                'margin-bottom': '5px',
                'margin-top': '35px',
                'color': colors['text']
            }
        ),

        # Nucleus
        html.Label('Nucleus'),
        html.Div(id='nucleus_widget_id', children=[
            dcc.Dropdown(
                id='nucleus_id'
            )],
            style={'margin-bottom': '10px', 'margin-top': '0px'}
        ),

        # Magnetic flux density
        html.Label(id='Magnetic_flux_density_output_container'),
        html.Div([dcc.Slider(
                id='magnetic_flux_density',
                min=0,
                max=30,
                step=0.1,
                value=9.4
            )],
            style={'margin-bottom': '10px', 'margin-top': '0px'}
        ),

        # Magic angle spinning
        html.Label(id='spinning_frequency_output_container'),
        html.Div(className='row', children=[
            dcc.Slider(
                className='col s6 m6 l6',
                id='spinning_frequency_in_kHz_coarse',
                min=0.0,
                max=110,
                step=5.0,
                value=0,
            ),
            dcc.Slider(
                className='col s6 m6 l6',
                id='spinning_frequency_in_kHz_fine',
                min=0,
                max=5,
                step=0.050,
                value=2.5,
            )],
            style={'margin-bottom': '10px', 'margin-top': '0px'}
        ),




        # dimension
        html.Div(className='row', children=[
            html.H6(className='col s6 m6 l6', children='Dimension parameters'),
            daq.BooleanSwitch(
                id-'ppm_switch',
                className='col s6 m6 l6',
                label='Show ppm',
                labelPosition='bottom',
                # size=40,
                style={
                    'margin-bottom': '0px',
                    'margin-top': '5px',
                    'color': colors['text']
                }
            )],
            style={
                'margin-bottom': '0px',
                'margin-top': '35px',
                'color': colors['text']
            }
        ),

        # Number of points
        html.Label(id='number_of_points_output_container'),
        html.Div([dcc.Slider(
                id='number_of_points',
                min=8,
                max=16,
                step=1,
                value=10,
                marks={8: '', 9: '', 10: '', 11: '', 12: '', 13: '', 14: '', 15: '', 16: ''}
            )],
            style={'margin-bottom': '10px', 'margin-top': '0px'}
        ),

        # Spectral width
        html.Label(id='frequency_bandwidth_output_container', style={'display': 'block'}),
        html.Div(className='row', children=[
            dcc.Slider(
                className='col s6 m6 l6',
                id='frequency_bandwidth_coarse',
                min=1.0,
                max=1000,
                step=10,
                value=100
            ),
            dcc.Slider(
                className='col s6 m6 l6',
                id='frequency_bandwidth_fine',
                min=0.0,
                max=10,
                step=0.050,
                value=5,
            )],
            style={'margin-bottom': '10px', 'margin-top': '0px'}
        ),

        # Reference offset
        html.Label(id='reference_offset_output_container'),
        html.Div(className='row', children=[
            dcc.Slider(
                className='col s6 m6 l6',
                id='reference_offset_coarse',
                min=0.0,
                max=100,
                step=5,
                value=0,
            ),
            dcc.Slider(
                className='col s6 m6 l6',
                id='reference_offset_fine',
                min=0,
                max=5,
                step=0.050,
                value=2.5,
            )],
        style={'margin-bottom': '10px', 'margin-top': '0px'}
        ),
    ]

