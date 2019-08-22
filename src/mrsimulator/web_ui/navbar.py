# -*- coding: utf-8 -*-
import dash_bootstrap_components as dbc
import dash_html_components as html
import base64


mrsimulator_logo = "docs/_static/mrsimulator-dark.png"
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
        dbc.Col(
            dbc.NavLink(
                "Documentation",
                href="https://mrsimulator.readthedocs.io/en/stable/",
                className="navbar-light expand-lg",
                # style={"color": "white", "hover": {"color": "#272727"}},
                # className="navbar navbar-expand-lg navbar-dark bg-dark",
            )
        ),
        dbc.Col(
            dbc.NavLink(
                "Github",
                href="https://github.com/DeepanshS/mrsimulator",
                # style={"color": "white"},
                # className="navbar navbar-expand-lg navbar-dark bg-dark",
            )
        ),
    ],
    className="ml-auto flex-nowrap mt-3 mt-md-0",
    align="center",
)

navbar = dbc.Navbar(
    [
        dbc.Row(
            [
                dbc.Col(
                    html.Img(
                        src="data:image/png;base64,{}".format(encoded_image.decode()),
                        height="40px",
                    )
                )
            ],
            align="center",
            no_gutters=True,
        ),
        dbc.NavbarToggler(id="navbar-toggler"),
        dbc.Collapse([nav_items], id="navbar-collapse", navbar=True),
    ],
    color="dark",
    sticky="top",
    dark=True,
    className="navbar navbar-expand-lg navbar-dark bg-dark",
)
