# -*- coding: utf-8 -*-
import dash_bootstrap_components as dbc
import dash_html_components as html
import base64


__title__ = "mrsimulator"
__sub_title__ = "A web application framework for NMR lineshape simulation."
bar = "docs/_static/mrsimulator-dark.png"
encoded_image = base64.b64encode(open(bar, "rb").read())

colors = {"background": "#e2e2e2", "text": "#085858"}
style_title = {"textAlign": "center", "background": "#e8e8e8", "text": "#085858"}


image_Card = dbc.CardImg(
    src="data:image/png;base64,{}".format(encoded_image.decode()),
    top=True,
    style={
        "background-color": "rgba(245, 245, 245, 0)",
        "opacity": "1.0",
        "border": "none",
    },
    className="card h-200",
)

titlebar = dbc.Jumbotron(
    [
        dbc.Container(
            [
                image_Card,
                html.Hr(className="my-2"),
                html.P(
                    children=__sub_title__,
                    style={"textAlign": "center", "color": colors["text"]},
                    className="lead",
                ),
            ],
            fluid=True,
        )
    ],
    fluid=True,
    className="display-1",
)
