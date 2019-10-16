# -*- coding: utf-8 -*-
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc
import dash_daq as daq

# tooltips provide additional information.
tooltips = []
tooltip_format = {"placement": "bottom", "delay": {"show": 250, "hide": 10}}
button_format = {"outline": True, "color": "dark", "size": "md"}

# info
info_button = dbc.Button(
    html.I(className="fas fa-info-circle"), id="indicator_status", **button_format
)
tooltips.append(dbc.Tooltip("Info", target="indicator_status", **tooltip_format))

# Advance settings
setting_button = dbc.Button(
    html.I(className="fas fa-cog"), id="advance-id", **button_format
)
tooltips.append(dbc.Tooltip("Advance setting", target="advance-id", **tooltip_format))

# modal page for advance setting
model_line_1 = dbc.Row(
    [
        dbc.Col(dbc.Label("Integration density")),
        dbc.Col(
            dbc.Input(
                type="number",
                value=70,
                # type="slider",
                min=0,
                max=4096,
                step=1,
                id="averaging_quality",
            )
        ),
    ]
)

model_line_1_info = dbc.Label(
    size="sm", id="number_of_averaging_points", style={"color": "#566573"}
)

model_line_2 = dbc.Row(
    [
        dbc.Col(dbc.Label("Integration volume")),
        dbc.Col(
            dcc.Dropdown(
                id="n_octants",
                options=[
                    {"label": "Octant", "value": 0},
                    {"label": "Hemisphere", "value": 1},
                    {"label": "Sphere", "value": 2},
                ],
                value=1,
                clearable=False,
            )
        ),
    ]
)

model_setting = dbc.Modal(
    [
        dbc.ModalHeader("Advance setting"),
        dbc.ModalBody(dbc.FormGroup([model_line_1, model_line_1_info, model_line_2])),
        dbc.ModalFooter(
            dbc.Button(
                "Close",
                id="close_setting",
                color="dark",
                className="ml-auto",
                outline=True,
            )
        ),
    ],
    id="modal_setting",
    role="document",
    # modalClassName="modal-dialog",
    className="modal-dialog",
)
# end of modal page for advance setting

# Show spectrum from individual isotopomers
decompose_button = dbc.Button(
    html.I(className="fas fa-align-left fas-shake"),
    id="decompose",
    **button_format,
    active=True,
)

tooltips.append(
    dbc.Tooltip(
        "Show simulation from individual isotopomers",
        target="decompose",
        placement="bottom",
        delay={"show": 250, "hide": 10},
    )
)

group_1_buttons = dbc.ButtonGroup(
    [decompose_button, info_button, setting_button, model_setting],
    className="btn-group mr-2",
)


# Download dataset
# As CSDM
csdm_button = dbc.Button(
    html.I(className="fas fa-download"), id="download_csdm_button", **button_format
)
tooltips.append(
    dbc.Tooltip(
        "Download file as .csdf", target="download_csdm_button", **tooltip_format
    )
)

group2_buttons = dbc.ButtonGroup([html.A(csdm_button, href="", id="download_csdm")])
toolbar = dbc.Col([group_1_buttons, group2_buttons, html.Div(tooltips)])
