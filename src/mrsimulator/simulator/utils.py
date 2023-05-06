import csdmpy as cp
import numpy as np

__author__ = "Deepansh Srivastava"


def fix_negative_axis(csdm_obj):
    reverse_index = [-i - 1 for i, x in enumerate(csdm_obj.x) if x.increment.value < 0]

    new_dvs = []
    for dv in csdm_obj.dependent_variables:
        y = dv.components
        y = y if reverse_index == [] else np.flip(y, axis=tuple(reverse_index))
        new_dv = cp.as_dependent_variable(array=y)
        new_dv.copy_metadata(dv)
        new_dvs.append(new_dv)

    new_dimensions = []
    for dim in csdm_obj.dimensions:
        if dim.increment.value < 0:
            coords = dim.coordinates
            array, unit = coords.value[::-1], coords.unit
            new_axis = cp.as_dimension(array=array, unit=unit)
            new_axis.copy_metadata(dim)
            np.allclose(coords.value[::-1], new_axis.coordinates.value)
            new_dimensions.append(new_axis)
        else:
            new_dimensions.append(dim)

    new_csdm = cp.CSDM(dimensions=new_dimensions, dependent_variables=new_dvs)
    new_csdm.copy_metadata(csdm_obj)

    return new_csdm
