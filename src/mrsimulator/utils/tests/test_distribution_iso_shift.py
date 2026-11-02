"""Tests for the isotropic chemical shift applied to tensor distribution line shapes
in the spectral fitting utilities, `_apply_iso_shift` and
`_generate_distribution_spectrum`.
"""
import csdmpy as cp
import mrsimulator.utils.spectral_fitting as sf
import numpy as np
import pytest
from mrsimulator.method.lib import BlochDecaySpectrum
from mrsimulator.models import CzjzekDistribution
from mrsimulator.models.utils import LineShapeKernel
from mrsimulator.simulator import ConfigSimulator


def _center_of_gravity(ppm, amp):
    return np.sum(ppm * amp) / np.sum(amp)


def _gaussian(ppm, center, sigma=8.0):
    return np.exp(-((ppm - center) ** 2) / (2 * sigma**2))


def _symmetric_ppm_axis(
    count, increment, complex_fft, origin_offset, reciprocal_offset
):
    """A LinearDimension covering a symmetric ppm range about zero."""
    offset = 0.0 if complex_fft else -increment * count / 2
    reciprocal = {}
    if reciprocal_offset:
        reciprocal = {"coordinates_offset": f"{-1 / (2 * increment)} s"}
    dim = cp.LinearDimension(
        count=count,
        increment=f"{increment} Hz",
        coordinates_offset=f"{offset} Hz",
        origin_offset=f"{origin_offset} Hz",
        complex_fft=complex_fft,
        reciprocal=reciprocal,
    )
    ppm_dim = dim.copy()
    ppm_dim.to("ppm", "nmr_frequency_ratio")
    return dim, ppm_dim.coordinates.value


AXIS_CONVENTIONS = [
    # (increment in Hz, complex_fft, reciprocal dimension carries an offset)
    (100.0, True, True),  # mrsimulator simulation output
    (100.0, True, False),  # hand-built csdmpy dimension, one-sided time axis
    (100.0, False, False),  # Bruker-like import, ascending axis
    (-100.0, True, True),
    (-100.0, False, False),  # Bruker-like import, descending axis
]


@pytest.mark.parametrize("increment, complex_fft, reciprocal_offset", AXIS_CONVENTIONS)
@pytest.mark.parametrize("larmor_sign", [-1.0, 1.0])
@pytest.mark.parametrize("iso_ppm", [60.0, 60.4])
def test_apply_iso_shift_direction(
    increment, complex_fft, reciprocal_offset, larmor_sign, iso_ppm
):
    """The FFT shift theorem must move a peak to +iso ppm regardless of
    - whether the dimension uses the complex_fft convention,
    - the sign of the frequency increment,
    - the sign of the Larmor frequency (gyromagnetic ratio) of the isotope,
    - whether the shift is an integer number of bins (60.0 vs 60.4 ppm at 1 ppm/bin),
    - whether the time axis returned by `fft()` is one-sided or centered.

    csdmpy's `CSDM.fft()` performs an inverse FFT when `complex_fft` is True and a
    forward FFT otherwise. The two use conjugate kernels, so a phase ramp of a fixed
    sign shifts the spectrum in opposite directions for the two conventions.
    """
    count, origin_offset = 1024, 100e6  # 1 ppm == 100 Hz
    # `Isotope.larmor_freq` returns -gamma * B0, so the sign depends on the isotope.
    larmor_freq = larmor_sign * origin_offset / 1e6

    dim, ppm = _symmetric_ppm_axis(
        count, increment, complex_fft, origin_offset, reciprocal_offset
    )
    amp = _gaussian(ppm, 0.0)
    obj = cp.CSDM(dimensions=[dim], dependent_variables=[cp.as_dependent_variable(amp)])
    assert abs(_center_of_gravity(ppm, amp)) < 1e-6

    shifted = sf._apply_iso_shift(
        obj.copy(), iso_shift_ppm=iso_ppm, larmor_freq=larmor_freq
    )
    y = shifted.y[0].components[0]
    msg = (
        f"increment={increment} Hz, complex_fft={complex_fft}, "
        f"reciprocal_offset={reciprocal_offset}, larmor_sign={larmor_sign}"
    )

    # the shifted spectrum must remain in pure absorption mode
    assert np.abs(y.imag).max() < 1e-6 * np.abs(y.real).max(), msg
    # and match the analytically shifted Gaussian
    assert np.allclose(y.real, _gaussian(ppm, iso_ppm), atol=1e-6), msg
    assert abs(_center_of_gravity(ppm, y.real) - iso_ppm) < 1e-3, msg


def _experiment_like(method, complex_fft, increment_sign):
    """Return an experiment CSDM object on the method's grid, either in mrsimulator's
    own convention (complex_fft=True) or in the convention of an imported dataset
    (complex_fft=False) with an ascending or descending frequency axis."""
    sd = method.spectral_dimensions[0]
    w_ref = method.channels[0].B0_to_ref_freq(method.magnetic_flux_density)
    if complex_fft:
        dim = sd.to_csdm_dimension()
        dim.origin_offset = f"{w_ref} Hz"
    else:
        increment = increment_sign * sd.spectral_width / sd.count
        dim = cp.LinearDimension(
            count=sd.count,
            increment=f"{increment} Hz",
            coordinates_offset=f"{-increment * sd.count / 2} Hz",
            origin_offset=f"{w_ref} Hz",
            complex_fft=False,
        )
    exp = cp.CSDM(
        dimensions=[dim],
        dependent_variables=[cp.as_dependent_variable(np.zeros(sd.count))],
    )
    exp.x[0].to("ppm", "nmr_frequency_ratio")
    return exp


@pytest.mark.parametrize("isotope", ["27Al", "29Si"])  # positive and negative gamma
@pytest.mark.parametrize(
    "complex_fft, increment_sign", [(True, 1), (False, 1), (False, -1)]
)
def test_generate_distribution_spectrum_iso_shift(isotope, complex_fft, increment_sign):
    """End-to-end: a Czjzek shielding distribution with a mean isotropic chemical shift
    of +60 ppm must appear at +60 ppm on the experiment's ppm axis, for a positive and
    a negative gyromagnetic ratio nucleus, and for experiments stored with or without
    the complex_fft convention. The spectral width is chosen so that the shift is not
    an integer number of bins."""
    iso_ppm = 60.0
    method = BlochDecaySpectrum(
        channels=[isotope],
        magnetic_flux_density=9.4,
        spectral_dimensions=[{"count": 1024, "spectral_width": 60e3}],
    )
    method.experiment = _experiment_like(method, complex_fft, increment_sign)

    pos = (np.linspace(-10, 10, 11), np.linspace(0, 1, 3))
    kernel = LineShapeKernel(
        method=method, pos=pos, config=ConfigSimulator(number_of_sidebands=1)
    )
    kernel.generate_lineshape(tensor_type="shielding")

    model = CzjzekDistribution(sigma=2.0, mean_isotropic_chemical_shift=iso_ppm)
    params = sf.make_LMFIT_params(spin_system_models=[model])
    spectrum = sf._generate_distribution_spectrum(params, kernel, [model])

    ppm = method.experiment.x[0].coordinates.value
    amp = spectrum.y[0].components[0].real
    assert amp.max() > 0
    # The shielding Czjzek distribution is symmetric in zeta, so the center of gravity
    # of the line shape is the mean isotropic chemical shift.
    assert (
        abs(_center_of_gravity(ppm, amp) - iso_ppm) < 0.1
    ), f"{isotope}, complex_fft={complex_fft}, increment_sign={increment_sign}"
