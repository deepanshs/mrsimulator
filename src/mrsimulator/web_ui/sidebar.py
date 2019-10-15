# -*- coding: utf-8 -*-
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc


colors = {"background": "#e2e2e2", "text": "#585858"}


filename_datetime = html.Div(
    [
        html.H5(id="filename_dataset"),
        html.H6(id="data_time", style={"textAlign": "left", "color": colors["text"]}),
    ]
)


upload_data = html.Div(
    [
        dcc.Upload(
            id="upload_data",
            children=html.Div(
                "Drag and drop, or click to select", style={"font-size": 16}
            ),
            style={
                "lineHeight": "60px",
                "borderWidth": "0.75px",
                "borderStyle": "dashed",
                "borderRadius": "7px",
                "textAlign": "center",
            },
            # Allow multiple files to be uploaded
            multiple=False,
            className="control-upload",
        ),
        dbc.Label(
            id="error_message",
            align="center",
            style={"font-size": 12, "align": "center"},
        ),
    ]
)


input_file = [html.Br(), upload_data, html.Hr(), filename_datetime]


sidebar = dbc.Card(
    dbc.CardBody(html.Div(input_file)), className="h-100 my-card", inverse=False
)
