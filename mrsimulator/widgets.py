

import dash_html_components as html
import dash_core_components as dcc


__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


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
                html.A(
                    id='download_link',
                    children='\u21E9 Download CSV'
                ),
                dcc.Graph(
                    id='nmr_spectrum',
                    figure={'data': []}
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
        html.Div([
            html.H6('Environment parameters')],
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
                marks={0: '+0 Hz', 50: '50 kHz', 110: '110 kHz'}
            ),
            dcc.Slider(
                className='col s6 m6 l6',
                id='spinning_frequency_in_kHz_fine',
                min=0,
                max=5,
                step=0.050,
                value=2.5,
                marks={2.5: '+2.5 kHz', 5: '+5 kHz'}
            )],
            style={'margin-bottom': '10px', 'margin-top': '0px'}
        ),




        # dimension
        html.Div(className='row', children=[
            html.H6(className='col s6 m6 l6', children='Dimension parameters'),
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
                marks={
                    8: '',
                    9: '',
                    10: '',
                    11: '',
                    12: '',
                    13: '',
                    14: '',
                    15: '',
                    16: ''
                }
            )],
            style={'margin-bottom': '10px', 'margin-top': '0px'}
        ),

        # Spectral width
        html.Label(id='frequency_bandwidth_output_container'),
        html.Div(className='row', children=[
            dcc.Slider(
                className='col s6 m6 l6',
                id='frequency_bandwidth_coarse',
                min=0.0,
                max=1000,
                step=50,
                value=100,
                marks={0: '0 Hz', 500: '0.5 MHz', 1000: '1 MHz'}
            ),
            dcc.Slider(
                className='col s6 m6 l6',
                id='frequency_bandwidth_fine',
                min=0.0,
                max=50,
                step=0.050,
                value=25,
                marks={0: '', 25: '+25 kHz', 50: '+50 kHz'}
            )],
            # style={'margin-bottom': '10px', 'margin-top': '0px'}
        ),

        # Reference offset
        html.Label(id='reference_offset_output_container'),
        html.Div(className='row', children=[
            dcc.Slider(
                className='col s6 m6 l6',
                id='reference_offset_coarse',
                min=0.0,
                max=100,
                step=10,
                value=0,
                marks={0: '0 Hz', 50: '50 kHz', 100: '100 kHz'},
            ),
            dcc.Slider(
                className='col s6 m6 l6',
                id='reference_offset_fine',
                min=0,
                max=10,
                step=0.050,
                value=0,
                marks={0: '', 5: '+5 kHz', 10: '+10 kHz'}
            )],
            # style={'margin-bottom': '10px', 'margin-top': '0px'}
        ),
    ]
