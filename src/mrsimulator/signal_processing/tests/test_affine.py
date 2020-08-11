# -*- coding: utf-8 -*-
import csdmpy as cp
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.affine_transformation as af
import numpy as np

# Creating a test CSDM object.
test_data = np.zeros((40, 40))
test_data[20, :] = 1
csdm_object = cp.CSDM(
    dependent_variables=[cp.as_dependent_variable(test_data)],
    dimensions=[
        cp.LinearDimension(count=40, increment="1 s", label="0", complex_fft=True),
        cp.LinearDimension(count=40, increment="1 K", label="1", complex_fft=True),
    ],
)

csdm_object2 = csdm_object.copy()


def test_shear_01():
    processor = sp.SignalProcessor(
        operations=[af.Shear(factor="-1 K/s", dim_index=1, normal=0)]
    )

    shear_data = processor.apply_operations(data=csdm_object)
    index = np.where(shear_data.dependent_variables[0].components[0] > 0.99999999)

    a = np.arange(40)
    assert np.allclose(index, [a, a])

    # complex_fft dim=0 to false
    csdm_object.dimensions[0].complex_fft = False
    csdm_object.dimensions[1].complex_fft = True
    shear_data = processor.apply_operations(data=csdm_object)
    index = np.where(shear_data.dependent_variables[0].components[0] > 0.99999999)

    a1 = np.arange(20)
    b1 = a1 + 20
    b = np.append(b1, a1)
    assert np.allclose(index, [a, b])

    # complex_fft dim=1 to false
    csdm_object.dimensions[0].complex_fft = True
    csdm_object.dimensions[1].complex_fft = False
    shear_data = processor.apply_operations(data=csdm_object)
    index = np.where(shear_data.dependent_variables[0].components[0] > 0.99999999)

    b = np.arange(40)
    b[1:] = a[::-1][:-1]
    assert np.allclose(index, [a, b])

    # both complex_fft set to false
    csdm_object.dimensions[0].complex_fft = False
    csdm_object.dimensions[1].complex_fft = False
    shear_data = processor.apply_operations(data=csdm_object)
    index = np.where(shear_data.dependent_variables[0].components[0] > 0.99999999)

    a1 = np.arange(21)[::-1]
    b1 = a1[1:-1] + 20
    b = np.append(a1, b1)
    assert np.allclose(index, [a, b])
