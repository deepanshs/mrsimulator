# -*- coding: utf-8 -*-
import dash_bootstrap_components as dbc
import dash_html_components as html
import base64
from os.path import join, split

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


folder = split(__file__)[0]
mrsimulator_logo = join(folder, "resource/mrsimulator-dark-featured.png")
encoded_image = base64.b64encode(open(mrsimulator_logo, "rb").read())

search_bar = dbc.Row(
    [
        dbc.Col(dbc.Input(type="search", placeholder="Search")),
        dbc.Col(dbc.Button("Search", color="primary", className="ml-2"), width="auto"),
    ],
    no_gutters=True,
    className="ml-auto flex-nowrap mt-3 mt-md-0",
    align="center",
)


nav_items = dbc.Row(
    [
        dbc.NavLink(
            "Documentation",
            href="https://mrsimulator.readthedocs.io/en/stable/",
            className="navbar-light expand-lg",
        ),
        dbc.NavLink(
            children=[html.I(className="fab fa-github-square"), "Github"],
            href="https://github.com/DeepanshS/mrsimulator",
        ),
    ],
    className="ml-auto",  # mt-3 mt-md-0 flex-nowrap
)

navbar_top = dbc.Navbar(
    [
        html.Div(
            [
                html.Img(
                    src="data:image/png;base64,{}".format(encoded_image.decode()),
                    height="70px",
                ),
                html.H5(
                    "A plotly-dash app",
                    style={"color": "white", "float": "center", "font-size": 16},
                ),
            ]
        ),
        dbc.NavbarToggler(id="navbar-toggler"),
        dbc.Collapse([nav_items], id="navbar-collapse", navbar=True),
    ],
    color="dark",
    sticky="top",
    fixed="top",
    dark=True,
    className="navbar navbar-expand-lg navbar-dark bg-dark",
)

navbar_bottom = dbc.Navbar(
    [
        dbc.Label("mrsimulator 2018-2019", style={"color": "white"}),
        html.I(className="fab fa-apple"),
    ],
    color="dark",
    sticky="bottom",
    # fixed="bottom",
    dark=True,
    # className="navbar navbar-expand-lg navbar-dark bg-dark",
)
