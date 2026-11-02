"""Apodization test"""
import csdmpy as cp
import numpy as np
from mrsimulator import signal_processor as sp
from mrsimulator import Simulator
from mrsimulator import SpinSystem
from mrsimulator.method.lib import BlochDecaySpectrum

from .test_signal_processor import setup_read_write

__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"

sim = Simulator()
the_site = {"isotope": "1H", "isotropic_chemical_shift": "0 ppm"}
the_spin_system = {"name": "site A", "sites": [the_site], "abundance": "80%"}
spin_system_object = SpinSystem.parse_dict_with_units(the_spin_system)
sim.spin_systems += [spin_system_object, spin_system_object, spin_system_object]
sim.config.decompose_spectrum = "spin_system"

sim.methods += [
    BlochDecaySpectrum(
        channels=["1H"],
        magnetic_flux_density=9.4,
        spectral_dimensions=[{"count": 65536, "spectral_width": 25000}],
    )
]
sim.run()
freqHz = sim.methods[0].spectral_dimensions[0].coordinates_Hz()


def test_scale():
    PS_0 = [sp.Scale(factor=10)]
    post_sim = sp.SignalProcessor(operations=PS_0)
    dataset = post_sim.apply_operations(dataset=sim.methods[0].simulation.copy())
    _, y0, y1, y2 = sim.methods[0].simulation.to_list()
    _, y0_, y1_, y2_ = dataset.to_list()
    # cast complex dataset
    assert np.allclose(y0_, y0 * 10), "Scaling failed"
    assert np.allclose(y1_, y1 * 10), "Scaling failed"
    assert np.allclose(y2_, y2 * 10), "Scaling failed"


def test_Lorentzian():
    PS_1 = [
        sp.IFFT(dim_index=0),
        sp.apodization.Exponential(FWHM="200 Hz", dim_index=0, dv_index=0),
        sp.FFT(dim_index=0),
    ]
    post_sim = sp.SignalProcessor(operations=PS_1)
    dataset = post_sim.apply_operations(dataset=sim.methods[0].simulation.copy())
    _, y0, y1, y2 = dataset.to_list()

    FWHM = 200
    test = (FWHM / 2) / (np.pi * (freqHz**2 + (FWHM / 2) ** 2))

    assert np.allclose(y1, y2)
    assert np.all(y0 != y1)
    assert np.allclose(
        test / test.sum(), y0 / y0.sum(), atol=1e-04
    ), "Lorentzian apodization amplitude failed"
    assert np.allclose(
        y1.sum(), y0.sum()
    ), "Area not conserved after Lorentzian apodization"


def test_Gaussian():
    FWHM = 200 * 2.354820045030949
    PS_2 = [
        sp.IFFT(dim_index=0),
        sp.apodization.Gaussian(FWHM=f"{FWHM} Hz", dim_index=0, dv_index=[0, 1]),
        sp.FFT(dim_index=0),
    ]

    PS_3 = [
        sp.IFFT(dim_index=0),
        sp.apodization.Gaussian(FWHM=f"{FWHM} Hz", dim_index=0, dv_index=None),
        sp.FFT(dim_index=0),
    ]

    post_sim = sp.SignalProcessor(operations=PS_2)
    dataset = post_sim.apply_operations(dataset=sim.methods[0].simulation.copy())
    _, y0, y1, _ = dataset.to_list()

    sigma = 200
    test = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((freqHz / sigma) ** 2) / 2)

    assert np.allclose(y0, y1), "Gaussian apodization on two dv are not equal."

    assert np.allclose(
        test / test.sum(), y0 / y0.sum(), atol=1e-04
    ), "Gaussian apodization amplitude failed"

    # test None for dv_index
    post_sim = sp.SignalProcessor(operations=PS_3)
    dataset = post_sim.apply_operations(dataset=sim.methods[0].simulation.copy())
    _, y0, y1, y2 = dataset.to_list()

    assert np.allclose(y0, y1), "Gaussian apodization on dv at 0 and 1 are unequal."
    assert np.allclose(y0, y2), "Gaussian apodization on dv at 0 and 2 are unequal."
    assert np.allclose(
        test / test.sum(), y0 / y0.sum(), atol=1e-04
    ), "Gaussian apodization amplitude failed"


def test_SkewedGaussian():
    # TODO: update this test for multiple skews and using np.convolve
    skew = 2
    FWHM = 200 * 2.354820045030949
    PS_2 = [
        sp.IFFT(dim_index=0),
        sp.apodization.SkewedGaussian(
            skew=skew, FWHM=f"{FWHM} Hz", dim_index=0, dv_index=[0, 1]
        ),
        sp.FFT(dim_index=0),
    ]

    post_sim = sp.SignalProcessor(operations=PS_2)
    dataset = post_sim.apply_operations(dataset=sim.methods[0].simulation.copy())
    _, y0, y1, _ = dataset.to_list()

    assert np.allclose(y0, y1), "Gaussian apodization on two dv are not equal."


def test_TopHat():
    test_dataset = cp.CSDM(
        dependent_variables=[cp.as_dependent_variable(np.ones(500))],
        dimensions=[cp.LinearDimension(500, "1 s")],
    )

    processor = sp.SignalProcessor()
    processor.operations = [
        sp.apodization.TopHat(rising_edge="100 s", falling_edge="400 s")
    ]

    rise_and_fall_dataset = processor.apply_operations(test_dataset.copy())
    rise_and_fall_should_be = np.zeros(500)
    rise_and_fall_should_be[100:400] = 1

    assert np.allclose(rise_and_fall_dataset.y[0].components, rise_and_fall_should_be)

    processor.operations = [sp.apodization.TopHat(rising_edge="100 s")]

    rise_only_dataset = processor.apply_operations(test_dataset.copy())
    rise_only_should_be = np.zeros(500)
    rise_only_should_be[100:] = 1

    assert np.allclose(rise_only_dataset.y[0].components, rise_only_should_be)

    processor.operations = [sp.apodization.TopHat(falling_edge="400 s")]

    fall_only_dataset = processor.apply_operations(test_dataset.copy())
    fall_only_should_be = np.zeros(500)
    fall_only_should_be[:400] = 1

    assert np.allclose(fall_only_dataset.y[0].components, fall_only_should_be)


def test_Mask():
    one_mask = np.zeros(shape=len(freqHz))

    PS_5 = [
        sp.IFFT(dim_index=0),
        sp.apodization.Mask(mask=one_mask, dim_index=0, dv_index=[0, 1]),
        sp.FFT(dim_index=0),
    ]

    post_sim = sp.SignalProcessor(operations=PS_5)
    dataset = post_sim.apply_operations(dataset=sim.methods[0].simulation.copy())
    _, y0, y1, _ = dataset.to_list()

    _, test_y0, test_y1, _ = sim.methods[0].simulation.to_list()

    nonzero_y0 = np.count_nonzero(y0)
    nonzero_y1 = np.count_nonzero(y1)

    assert np.allclose(y0, y1), "Mask on two dv are not equal."

    assert np.allclose(
        nonzero_y0, nonzero_y1, atol=1e-04
    ), "Mask apodization amplitude failed"


def test_scale_class():
    # direct initialization
    a = sp.Scale(factor=200)

    assert a.factor == 200
    assert a.property_units == {}

    # class to dict with units
    dict_ = a.json()

    assert dict_ == {
        "function": "Scale",
        "factor": 200.0,
    }

    # read from dictionary
    b = sp.Scale.parse_dict_with_units(dict_)

    assert a == b


def test_TopHat_class():
    a = sp.apodization.TopHat(rising_edge="1 s")
    assert a.property_units == {"rising_edge": "s"}

    a = sp.apodization.TopHat(falling_edge="1 s")
    assert a.property_units == {"falling_edge": "s"}


def test_Exponential_class():
    # direct initialization
    a = sp.apodization.Exponential(FWHM="200 s", dim_index=0, dv_index=0)

    assert a.FWHM == 200
    assert a.property_units == {"FWHM": "s"}
    assert a.dim_index == 0
    assert a.dv_index == 0

    py_dict = {
        "function": "apodization",
        "type": "Exponential",
        "FWHM": "200.0 s",
        "dim_index": 0,
        "dv_index": 0,
    }
    setup_read_write(a, py_dict, sp.apodization.Exponential)


def test_Gaussian_class():
    # direct initialization
    a = sp.apodization.Gaussian(FWHM="200 km/s", dim_index=0, dv_index=0)

    assert a.FWHM == 200
    assert a.property_units == {"FWHM": "km / s"}
    assert a.dim_index == 0
    assert a.dv_index == 0

    py_dict = {
        "function": "apodization",
        "type": "Gaussian",
        "FWHM": "200.0 km / s",
        "dim_index": 0,
        "dv_index": 0,
    }
    setup_read_write(a, py_dict, sp.apodization.Gaussian)


def test_2D_area():
    data = np.zeros((256, 128), dtype=float)
    data[128, 64] = 1.0
    csdm_obj = cp.as_csdm(data)
    csdm_obj.x[0].period = "256"
    csdm_obj.x[1].period = "128"

    # test00
    PS = [
        sp.IFFT(dim_index=(0, 1)),
        sp.apodization.Gaussian(FWHM="35", dim_index=0),
        sp.apodization.Gaussian(FWHM="55", dim_index=1),
        sp.FFT(dim_index=(0, 1)),
    ]

    post_sim = sp.SignalProcessor(operations=PS)
    dataset_new = post_sim.apply_operations(dataset=csdm_obj.copy())
    _, __, y1 = dataset_new.to_list()

    assert np.allclose(y1.sum(), data.sum())

    # test01
    PS = [
        sp.IFFT(dim_index=(0, 1)),
        sp.apodization.Gaussian(FWHM="35", dim_index=0),
        sp.apodization.Exponential(FWHM="55", dim_index=1),
        sp.FFT(dim_index=(0, 1)),
    ]

    post_sim = sp.SignalProcessor(operations=PS)
    dataset_new = post_sim.apply_operations(dataset=csdm_obj.copy())
    _, __, y1 = dataset_new.to_list()

    assert np.allclose(y1.sum(), data.sum())


def test_MultiDimensionApodization_class():
    a = sp.apodization.MultiDimensionApodization()
    assert a.function == "apodization"


def test_TopHat_mixed_units():
    """Each edge keeps its own unit, and json keeps both units."""
    test_dataset = cp.CSDM(
        dependent_variables=[cp.as_dependent_variable(np.ones(500))],
        dimensions=[cp.LinearDimension(500, "1 ms")],
    )
    top_hat = sp.apodization.TopHat(rising_edge="0.1 s", falling_edge="400 ms")
    assert str(top_hat.property_units["rising_edge"]) == "s"
    assert str(top_hat.property_units["falling_edge"]) == "ms"
    assert top_hat.json()["rising_edge"] == "0.1 s"
    assert top_hat.json()["falling_edge"] == "400.0 ms"

    processor = sp.SignalProcessor(operations=[top_hat])
    out = processor.apply_operations(test_dataset.copy())
    expected = np.zeros(500)
    expected[100:400] = 1
    assert np.allclose(out.y[0].components[0], expected)

    # no edges selects the full window
    processor = sp.SignalProcessor(operations=[sp.apodization.TopHat()])
    out = processor.apply_operations(test_dataset.copy())
    assert np.allclose(out.y[0].components[0], 1)


def _time_domain_fid(n=2**16, dt_us=10.0):
    """Two-sided constant signal, whose spectrum is the Fourier transform of the
    apodization window."""
    return cp.CSDM(
        dimensions=[cp.LinearDimension(n, f"{dt_us} us", complex_fft=True)],
        dependent_variables=[cp.as_dependent_variable(np.ones(n, dtype=complex))],
    )


def _spectrum(op):
    """Spectrum of the simulated delta line (1H at 0 ppm) after apodization in the
    time domain; (frequency in Hz, complex amplitude)."""
    processor = sp.SignalProcessor(operations=[sp.IFFT(), op, sp.FFT()])
    out = processor.apply_operations(dataset=sim.methods[0].simulation.copy())
    return freqHz, out.y[0].components[0]


def _fwhm(x, y):
    above = x[y >= y.max() / 2]
    return above.max() - above.min() + abs(x[1] - x[0])


def test_SkewedGaussian_zero_skew_is_Gaussian():
    """With zero skew, the skewed Gaussian window is the Gaussian window of the same
    FWHM, with unit amplitude at time zero."""
    fid = _time_domain_fid()
    gauss = sp.apodization.Gaussian(FWHM="100 Hz").operate(fid.copy())
    skewed = sp.apodization.SkewedGaussian(FWHM="100 Hz", skew=0).operate(fid.copy())
    assert np.allclose(skewed.y[0].components[0], gauss.y[0].components[0])


def test_SkewedGaussian_unit_independent():
    """The window does not depend on the unit the FWHM is given in."""
    fid = _time_domain_fid()
    a = sp.apodization.SkewedGaussian(FWHM="100 Hz", skew=2).operate(fid.copy())
    b = sp.apodization.SkewedGaussian(FWHM="0.1 kHz", skew=2).operate(fid.copy())
    assert np.allclose(a.y[0].components[0], b.y[0].components[0])


def test_SkewedGaussian_lineshape():
    """The apodized delta line is a skew normal line: real, of area equal to the
    Gaussian line, with the FWHM set by the FWHM argument at zero skew, and skewed
    towards higher frequencies for a positive skew."""
    x, y_gauss = _spectrum(sp.apodization.Gaussian(FWHM="100 Hz"))
    x, y0 = _spectrum(sp.apodization.SkewedGaussian(FWHM="100 Hz", skew=0))
    x, y3 = _spectrum(sp.apodization.SkewedGaussian(FWHM="100 Hz", skew=3))
    x, ym3 = _spectrum(sp.apodization.SkewedGaussian(FWHM="100 Hz", skew=-3))

    assert np.allclose(_fwhm(x, y0.real), 100.0, rtol=0.03)
    for y in [y0, y3, ym3]:
        assert np.allclose(y.imag, 0, atol=1e-6 * y.real.max())
        assert np.all(y.real > -1e-6 * y.real.max())
        assert np.allclose(y.real.sum(), y_gauss.real.sum())

    def third_moment(y):
        p = y.real / y.real.sum()
        mean = (p * x).sum()
        return (p * (x - mean) ** 3).sum()

    assert third_moment(y3) > 0
    assert third_moment(ym3) < 0
    assert np.allclose(third_moment(y3), -third_moment(ym3))
    assert abs(third_moment(y0)) < 1e-3 * abs(third_moment(y3))


def test_SkewedGaussian_real_dataset():
    """A complex window applied to a real-valued dataset returns the complex product
    instead of failing on the in-place cast."""
    n = 256
    dataset = cp.CSDM(
        dimensions=[cp.LinearDimension(n, "10 us", complex_fft=True)],
        dependent_variables=[cp.as_dependent_variable(np.ones(n))],
    )
    assert dataset.y[0].components.dtype == np.float64
    op = sp.apodization.SkewedGaussian(FWHM="100 Hz", skew=2)
    out = op.operate(dataset.copy())
    expected = op.fn(dataset.x[0].coordinates)
    assert np.iscomplexobj(out.y[0].components)
    assert np.allclose(out.y[0].components[0], expected)

    # a real window keeps the real dtype
    out = sp.apodization.Gaussian(FWHM="100 Hz").operate(dataset.copy())
    assert out.y[0].components.dtype == np.float64
