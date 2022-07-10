__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


def get_spectral_dimensions(csdm_object, units=False):
    """Extract the count, spectral_width, and reference_offset parameters, associated
    with the spectral dimensions of the method, from the CSDM dimension objects.

    Args:
        csdm_object: A CSDM object holding the measurement dataset.

    Returns:
        A list of dict objects, where each dict contains the count, spectral_width, and
        reference_offset.
    """
    result = []
    for dim in csdm_object.dimensions:
        count = dim.count
        increment = dim.increment.to("Hz").value
        ref = dim.coordinates_offset.to("Hz").value
        sw = count * increment

        even = count % 2 == 0
        complex_co = ref + sw / 2.0
        complex_co -= 0 if even else increment / 2.0
        co = ref if dim.complex_fft else complex_co

        # if sw < 0:
        #     sw = -sw
        #     co += -increment if even else 0

        dim_i = {}
        dim_i["count"] = dim.count
        dim_i["spectral_width"] = sw if not units else f"{sw} Hz"
        dim_i["reference_offset"] = co if not units else f"{co} Hz"

        if dim.label != "":
            dim_i["label"] = dim.label

        if dim.origin_offset not in ["", None, 0]:
            oo = dim.origin_offset.to("Hz").value
            dim_i["origin_offset"] = oo if not units else f"{oo} Hz"
        result.append(dim_i)

    return result[::-1]


def flatten_dict(obj, previous_key=None):
    """Flatten a nested dictionary with keys as obj1.obj2... and so on"""
    result = {}
    for k, v in obj.items():
        if not isinstance(v, dict):
            key = f"{previous_key}.{k}" if previous_key is not None else k
            result.update({key: v})
        else:
            result.update(**flatten_dict(v, previous_key=k))

    return result


# def plotly_scatter_obj(data: cp.CSDM = None, label: str = None):
#     """Plotly figure dict.
#     Example:

#         # >>> fig_obj = plotly_scatter_obj(experiment, label="experiment")
#         # >>> fig_obj["data"] += plotly_scatter_obj(best_fit)["data"]
#         # >>> fig_obj["data"] += plotly_scatter_obj(residuals, label="res")["data"]
#         # >>> fig = go.Figure(fig_obj)
#         # >>> fig
#     """
#     data_ = [
#         go.Scatter(
#             x=data.x[0].coordinates.value,
#             y=item.components[0],
#             mode="lines",
#             opacity=0.75,
#             name=item.name if label is None else label,
#         )
#         for item in data.y
#     ]
#     layout_ = go.Layout(
#         xaxis=dict(
#             title=data.x[0].axis_label,
#             ticks="outside",
#             showline=True,
#             autorange="reversed",
#             zeroline=False,
#         ),
#         yaxis=dict(
#             ticks="outside",
#             showline=True,
#             zeroline=False,
#             autorange=True,
#         ),
#         margin={"l": 60, "b": 45, "t": 50, "r": 50},
#         width=800,
#         height=400,
#     )

#     return {"data": data_, "layout": layout_}
