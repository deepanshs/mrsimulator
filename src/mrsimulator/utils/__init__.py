# -*- coding: utf-8 -*-

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


def get_spectral_dimensions(csdm_object):
    """
    Extract the count, spectral_width, and reference_offset parameters, associated with
    the spectral dimensions of the method, from the CSDM dimension objects.

    Args:
        csdm_object: A CSDM object holding the measurement dataset.

    Returns:
        A list of dict objects, where each dict containts the count, spectral_width, and
        reference_offset.
    """
    result = []
    for dim in csdm_object.dimensions:
        count = dim.count
        increment = dim.increment.to("Hz").value
        coordinates_offset = dim.coordinates_offset.to("Hz").value
        spectral_width = count * increment

        dim_i = {}
        dim_i["count"] = dim.count
        dim_i["spectral_width"] = spectral_width
        if dim.complex_fft is True:
            dim_i["reference_offset"] = coordinates_offset
        else:
            dim_i["reference_offset"] = coordinates_offset + spectral_width / 2.0

        if dim.label != "":
            dim_i["label"] = dim.label
        result.append(dim_i)

    return result[::-1]
