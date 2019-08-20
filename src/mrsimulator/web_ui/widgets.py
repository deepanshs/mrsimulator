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
        html.Br(),
        dbc.Row(
            [
                dbc.Col(spectrum_body, xs=12, sm=12, md=12, lg=7, xl=7),
                dbc.Col(dimension_body),
            ],
            no_gutters=False,
        ),
        html.Br(),
        # dbc.Jumbotron(),
        html.Br(),
    ],
    fluid=True,
    # className="mt-4",
)

input_file = (
    html.Div(
        className="row",
        children=[
            html.Div(
                className="col s12 m7 l7",
                children=[
                    html.H5(id="filename_dataset"),
                    html.H6(
                        id="data_time",
                        style={"textAlign": "left", "color": colors["text"]},
                    ),
                ],
            ),
            html.Div(
                className="col s12 m5 l5",
                children=[
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
                    html.Label(
                        id="error_message",
                        # 'color': 'red'}
                        style={"textAlign": "center"},
                    ),
                ],
            ),
        ],
    ),
)


def table_css():
    return {
        "selector": ".dash-cell div.dash-cell-value",
        "rule": {
            "display": "inline",
            "white-space": "inherit",
            "overflow": "inherit",
            "text-overflow": "inherit",
            "overflowX": "scroll",
            "border": "thin lightgrey solid",
        },
    }


def display_isotopomers(isotope, isotopomers_object):
    columns_ = ["index", "site", "isotropic chemical shift", "anisotropy", "asymmetry"]

    child = []
    for j, isotopomer in enumerate(isotopomers_object):
        sites_ = []
        for i, site in enumerate(isotopomer.sites):
            if site.isotope == isotope:
                site_ = {}
                site_["index"] = i
                # site_["isotope"] = site.isotope
                # site_["isotropic chemical shift"] = site.isotropic_chemical_shift
                # site_["anisotropy"] = site.shielding_symmetric.anisotropy
                # site_["asymmetry"] = site.shielding_symmetric.asymmetry

                # ravel_.append(site['isotope'])
                # ravel_.append(site['isotropic_chemical_shift'])
                # ravel_.append(site['shielding_symmetric']['anisotropy'])
                # ravel_.append(site['shielding_symmetric']['asymmetry'])
                sites_.append(site_)

                # print(sites_)

                # child.append(
                #     html.Div(
                #         className="row",
                #         children=[
                #             html.Div(
                #                 className="card-panel hoverable col s12 m6 l6",
                #                 children=[
                #                     html.H6(
                #                         children="".join(["Isotopomer ", str(j)]),
                #                         style={
                #                             "textAlign": "left",
                #                             "color": colors["text"],
                #                         },
                #                     ),
                #                     dash_table.DataTable(
                #                         style_data={"whiteSpace": "normal"},
                #                         css=[table_css()],
                #                         style_as_list_view=True,
                #                         style_cell={
                #                             "textAlign": "left",
                #                             "padding": "5px",
                #                         },
                #                         style_header={
                #                             "backgroundColor": "white",
                #                             "fontWeight": "bold",
                #                             "fontSize": 12,
                #                             "color": "#585858",
                #                         },
                #                         style_cell_conditional=[
                #                             {
                #                                 "if": {"row_index": "odd"},
                #                                 "backgroundColor": "#f8f8f8",
                #                             }
                #                         ],
                #                         columns=[
                #                             {"name": i, "id": i} for i in columns_
                #                         ],
                #                         data=sites_,
                #                     ),
                #                 ],
                #             )
                #         ],
                #     )
                # )

    # abundance = isotopomer['abundance']
    return child
