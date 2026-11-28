import csdmpy as cp
import numpy as np
from mrsimulator import signal_processor as sp

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


# Creating a test CSDM object.
test_data = np.zeros((40, 40))
test_data[20, :] = 1
csdm_object = cp.CSDM(
    dependent_variables=[cp.as_dependent_variable(test_data)],
    dimensions=[
        cp.LinearDimension(count=40, increment="1 s", label="0", complex_fft=True),
        cp.Dimension(
            type="linear", count=40, increment="1 K", label="1", complex_fft=True
        ),
    ],
)

csdm_object2 = csdm_object.copy()


def test_shear_01():
    processor = sp.SignalProcessor(
        operations=[
            sp.IFFT(dim_index=1),
            sp.affine.Shear(factor="-1 K/s", dim_index=1, parallel=0),
            sp.FFT(dim_index=1),
        ]
    )

    shear_dataset = processor.apply_operations(dataset=csdm_object)
    index = np.where(shear_dataset.y[0].components[0] > 0.99999999)

    a = np.arange(40)
    assert np.allclose(index, [a, a])

    # complex_fft dim=0 to false
    csdm_object.x[0].complex_fft = False
    csdm_object.x[1].complex_fft = True
    shear_dataset = processor.apply_operations(dataset=csdm_object)
    index = np.where(shear_dataset.y[0].components[0] > 0.99999999)

    a1 = np.arange(20)
    b1 = a1 + 20
    b = np.append(b1, a1)
    assert np.allclose(index, [a, b])

    # complex_fft dim=1 to false
    csdm_object.x[0].complex_fft = True
    csdm_object.x[1].complex_fft = False
    shear_dataset = processor.apply_operations(dataset=csdm_object)
    index = np.where(shear_dataset.y[0].components[0] > 0.99999999)

    b = np.arange(40)
    b[1:] = a[::-1][:-1]
    assert np.allclose(index, [a, b])

    # both complex_fft set to false
    csdm_object.x[0].complex_fft = False
    csdm_object.x[1].complex_fft = False
    shear_dataset = processor.apply_operations(dataset=csdm_object)
    index = np.where(shear_dataset.y[0].components[0] > 0.99999999)

    a1 = np.arange(21)[::-1]
    b1 = a1[1:-1] + 20
    b = np.append(a1, b1)
    assert np.allclose(index, [a, b])


def test_serialization_and_parse():
    processor = sp.SignalProcessor(
        operations=[
            sp.IFFT(dim_index=1),
            sp.affine.Shear(factor="-1 K/s", dim_index=1, parallel=0),
            sp.FFT(dim_index=1),
        ]
    )
    serialize = processor.json()

    expected = {
        "operations": [
            {"dim_index": 1, "function": "IFFT"},
            {
                "dim_index": 1,
                "factor": "-1.0 K / s",
                "function": "affine",
                "parallel": 0,
                "type": "Shear",
            },
            {"dim_index": 1, "function": "FFT"},
        ],
    }

    assert serialize == expected

    recovered = sp.SignalProcessor.parse_dict_with_units(serialize)
    assert recovered == processor


def test_scale():
    c_inc = csdm_object.x[1].increment.value
    c_off = csdm_object.x[1].coordinates_offset.value

    processor = sp.SignalProcessor(operations=[sp.affine.Scale(factor=2, dim_index=1)])
    scaled_dataset = processor.apply_operations(dataset=csdm_object)

    s_inc = scaled_dataset.x[1].increment.value
    s_off = scaled_dataset.x[1].coordinates_offset.value

    assert s_inc == 2 * c_inc
    assert s_off == 2 * c_off


def _shear_phase(unit_x, unit_y, factor, dtype=complex):
    """Phase / (-2 pi) at x = 1, y = 1, in the units of the dimensions."""
    dataset = cp.CSDM(
        dimensions=[
            cp.LinearDimension(count=8, increment=f"1 {unit_x}"),
            cp.LinearDimension(count=4, increment=f"1 {unit_y}"),
        ],
        dependent_variables=[cp.as_dependent_variable(np.ones((4, 8), dtype=dtype))],
    )
    out = sp.affine.Shear(factor=factor, dim_index=0, parallel=1).operate(dataset)
    return np.angle(out.y[0].components[0][1, 1]) / (-2 * np.pi)


def test_shear_units_of_dimensions():
    """The shear phase is factor * x * y, with x * y in the reciprocal of the unit
    of the factor, for dimensions in any unit."""
    assert np.allclose(_shear_phase("s", "Hz", 0.01), 0.01)
    assert np.allclose(_shear_phase("ms", "Hz", 0.01), 1e-5)
    assert np.allclose(_shear_phase("s", "kHz", 1e-5), 0.01)
    assert np.allclose(_shear_phase("ms", "kHz", 0.01), 0.01)
    # factor of 0.01 ms/s = 1e-5
    assert np.allclose(_shear_phase("s", "kHz", "0.01 ms/s"), 0.01)


def test_shear_real_dataset():
    """The shear of a real-valued dataset is complex."""
    assert np.allclose(_shear_phase("s", "Hz", 0.01, dtype=float), 0.01)
