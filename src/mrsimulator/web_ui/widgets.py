# -*- coding: utf-8 -*-
import dash_core_components as dcc
import dash_html_components as html
from mrsimulator.web_ui import dimension
from mrsimulator.web_ui.dimension import dimension_body
from mrsimulator.web_ui.graph import spectrum_body
import dash_table
import dash_bootstrap_components as dbc


__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


colors = {"background": "#e2e2e2", "text": "#585858"}

main_body = dbc.Row(
    [
        dbc.Col(spectrum_body, xs=12, sm=12, md=7, lg=7, xl=7),
        dbc.Col(dimension_body, xs=12, sm=12, md=5, lg=5, xl=5),
    ],
    no_gutters=False,
)


def get_isotopomers(isotope, isotopomers_object):
    list_isotopomers = [
        dbc.ListGroup(
            [
                dbc.ListGroupItem(
                    [
                        dbc.ListGroupItemHeading(f"Isotopomer {i+1}"),
                        dbc.ListGroupItem(
                            f"{j+1} {site.isotope} Site", n_clicks=0, action=True
                        ),
                    ]
                )
                for i, isotopomer in enumerate(isotopomers_object)
                for j, site in enumerate(isotopomer.sites)
            ]
        )
    ]
    return html.Div(list_isotopomers)
