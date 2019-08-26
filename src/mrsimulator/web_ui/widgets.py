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

main_body = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(spectrum_body, xs=12, sm=12, md=12, lg=7, xl=7),
                dbc.Col(dimension_body),
            ],
            no_gutters=False,
        )
    ],
    # fluid=True,
    # className="w3-container w3-teal",
)

filename_datetime = dbc.Col(
    [
        html.H5(id="filename_dataset"),
        html.P(id="data_time", style={"textAlign": "left", "color": colors["text"]}),
    ]
)

upload_data = dbc.Col(
    [
        dcc.Upload(
            id="upload_data",
            children=html.Div(["Drag and Drop or ", html.A("Select File")]),
            style={
                "width": "100%",
                "height": "50px",
                "lineHeight": "50px",
                "borderWidth": "1px",
                "borderStyle": "dashed",
                "borderRadius": "5px",
                "textAlign": "center",
                "margin": "1px",
            },
            # Allow multiple files to be uploaded
            multiple=False,
        ),
        dbc.Label(id="error_message", align="center"),
    ]
)

input_file = dbc.Container(dbc.Row([filename_datetime, upload_data]), fluid=True)


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


# def get_isotopomers(isotope, isotopomers_object):

#     cards = []
#     total_sites = 0
#     for j, isotopomer in enumerate(isotopomers_object):
#         number_of_sites = len(isotopomer.sites)
#         total_sites += number_of_sites
#         for i, site in enumerate(isotopomer.sites):
#             if site.isotope == isotope:
#                 card_body = [
#                     html.H6(f"{isotope} site index {i}", className="card-title"),
#                     dbc.Row(
#                         [
#                             dbc.Col(dbc.Label(
#                                  f"Isotropic chemical shift", size="sm"
#                             )),
#                             dbc.Col(
#                                 dbc.Label(
#                                     "{0} {1}".format(
#                                         site.isotropic_chemical_shift,
#                                         site.property_units["isotropic_chemical_shift"],
#                                     ),
#                                     className="card-text",
#                                     size="sm",
#                                     align="right",
#                                 )
#                             ),
#                         ],
#                         no_gutters=True,
#                     ),
#                     dbc.Row(
#                         [
#                             dbc.Col(dbc.Label(f"Shielding anisotropy", size="sm")),
#                             dbc.Col(
#                                 dbc.Label(
#                                     "{0} {1}".format(
#                                         site.shielding_symmetric.anisotropy,
#                                         site.shielding_symmetric.property_units[
#                                             "anisotropy"
#                                         ],
#                                     ),
#                                     className="card-text",
#                                     size="sm",
#                                     align="right",
#                                 )
#                             ),
#                         ],
#                         no_gutters=True,
#                     ),
#                     dbc.Row(
#                         [
#                             dbc.Col(dbc.Label(f"Shielding asymmetry", size="sm")),
#                             dbc.Col(
#                                 dbc.Label(
#                                     f"{site.shielding_symmetric.asymmetry}",
#                                     className="card-text",
#                                     size="sm",
#                                     align="right",
#                                 )
#                             ),
#                         ],
#                         no_gutters=True,
#                     ),
#                 ]
#                 card = [
#                     dbc.CardHeader(f"Isotopomer index {j}"),
#                     dbc.CardBody(card_body),
#                 ]
#                 cards.append(dbc.Card(card))

#     cards_deck = dbc.Jumbotron(dbc.CardDeck(cards))
#     return cards_deck, total_sites
