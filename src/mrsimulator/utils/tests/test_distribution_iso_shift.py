"""Tests for the isotropic chemical shift applied to tensor distribution line shapes
in the spectral fitting utilities, `_apply_iso_shift` and
`_generate_distribution_spectrum`.
"""
import csdmpy as cp
import mrsimulator.utils.spectral_fitting as sf
import numpy as np
import pytest
from mrsimulator.method.lib import BlochDecayCTSpectrum
from mrsimulator.method.lib import BlochDecaySpectrum
from mrsimulator.models import CzjzekDistribution
from mrsimulator.models.utils import LineShapeKernel
from mrsimulator.simulator import ConfigSimulator
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.utils import get_spectral_dimensions


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


def _experiment(isotope, B0, count, spectral_width, complex_fft, increment_sign):
    """Return an experiment CSDM object (in Hz), either in mrsimulator's own convention
    (complex_fft=True) or in the convention of an imported dataset (complex_fft=False)
    with an ascending or descending frequency axis."""
    w_ref = Isotope(symbol=isotope).B0_to_ref_freq(B0)
    increment = increment_sign * spectral_width / count
    offset = 0.0 if complex_fft else -increment * count / 2
    dim = cp.LinearDimension(
        count=count,
        increment=f"{increment} Hz",
        coordinates_offset=f"{offset} Hz",
        origin_offset=f"{w_ref} Hz",
        complex_fft=complex_fft,
    )
    return cp.CSDM(
        dimensions=[dim],
        dependent_variables=[cp.as_dependent_variable(np.zeros(count))],
    )


def _method_for_experiment(method_class, isotope, exp, **kwargs):
    """Documented fitting flow: the method grid is taken from the experiment."""
    method = method_class(
        channels=[isotope],
        magnetic_flux_density=9.4,
        spectral_dimensions=get_spectral_dimensions(exp),
        **kwargs,
    )
    exp.x[0].to("ppm", "nmr_frequency_ratio")
    method.experiment = exp
    return method


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
    exp = _experiment(isotope, 9.4, 1024, 60e3, complex_fft, increment_sign)
    method = _method_for_experiment(BlochDecaySpectrum, isotope, exp)

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


def _quad_distribution_spectrum(complex_fft, increment_sign):
    """Quadrupolar Czjzek distribution spectrum of 27Al (central transition, static)
    with a mean isotropic chemical shift of 60 ppm. Returns (ppm, amplitude)."""
    exp = _experiment("27Al", 9.4, 1024, 100e3, complex_fft, increment_sign)
    method = _method_for_experiment(
        BlochDecayCTSpectrum, "27Al", exp, rotor_frequency=0
    )
    pos = (np.linspace(0, 8e6, 17), np.linspace(0, 1, 6))
    kernel = LineShapeKernel(
        method=method, pos=pos, config=ConfigSimulator(number_of_sidebands=1)
    )
    kernel.generate_lineshape(tensor_type="quadrupolar")

    model = CzjzekDistribution(sigma=1.5e6, mean_isotropic_chemical_shift=60.0)
    params = sf.make_LMFIT_params(spin_system_models=[model])
    spectrum = sf._generate_distribution_spectrum(params, kernel, [model])
    ppm = method.experiment.x[0].coordinates.value
    return ppm, spectrum.y[0].components[0].real


def test_quadrupolar_distribution_spectrum_axis_conventions():
    """An asymmetric (second-order quadrupolar) distribution line shape must be
    identical on the ppm axis for experiments stored with or without the complex_fft
    convention, and with an ascending or descending frequency axis."""
    ppm_ref, amp_ref = _quad_distribution_spectrum(True, 1)
    amp_ref /= amp_ref.max()

    # the line shape is asymmetric about its maximum, so a mirrored kernel is detected
    peak = ppm_ref[np.argmax(amp_ref)]
    mirrored = np.interp(2 * peak - ppm_ref, ppm_ref, amp_ref)
    assert np.abs(mirrored - amp_ref).max() > 0.2

    for complex_fft, increment_sign in [(False, 1), (False, -1)]:
        ppm, amp = _quad_distribution_spectrum(complex_fft, increment_sign)
        order = np.argsort(ppm)
        amp_i = np.interp(ppm_ref, ppm[order], amp[order] / amp.max())
        assert np.abs(amp_i - amp_ref).max() < 1e-2, (complex_fft, increment_sign)


def test_generate_distribution_spectrum_mismatched_grid():
    """A method whose spectral grid does not match the experiment grid (here, an
    ascending method grid with a descending experiment axis) must raise an error
    instead of silently producing a mirrored spectrum."""
    exp = _experiment("27Al", 9.4, 1024, 100e3, False, -1)
    exp.x[0].to("ppm", "nmr_frequency_ratio")
    method = BlochDecayCTSpectrum(
        channels=["27Al"],
        magnetic_flux_density=9.4,
        rotor_frequency=0,
        spectral_dimensions=[{"count": 1024, "spectral_width": 100e3}],
        experiment=exp,
    )
    pos = (np.linspace(0, 8e6, 5), np.linspace(0, 1, 3))
    kernel = LineShapeKernel(
        method=method, pos=pos, config=ConfigSimulator(number_of_sidebands=1)
    )
    kernel.generate_lineshape(tensor_type="quadrupolar")

    model = CzjzekDistribution(sigma=1.5e6, mean_isotropic_chemical_shift=60.0)
    params = sf.make_LMFIT_params(spin_system_models=[model])
    with pytest.raises(ValueError, match="get_spectral_dimensions"):
        sf._generate_distribution_spectrum(params, kernel, [model])


def test_kernel_grid_checked_once_per_experiment(monkeypatch):
    """The kernel grid check runs once per (kernel, experiment) pair, not on every
    evaluation of the distribution spectrum, and again when the experiment changes."""
    exp = _experiment("27Al", 9.4, 256, 100e3, True, 1)
    method = _method_for_experiment(
        BlochDecayCTSpectrum, "27Al", exp, rotor_frequency=0
    )
    pos = (np.linspace(0, 8e6, 5), np.linspace(0, 1, 3))
    kernel = LineShapeKernel(
        method=method, pos=pos, config=ConfigSimulator(number_of_sidebands=1)
    )
    kernel.generate_lineshape(tensor_type="quadrupolar")
    model = CzjzekDistribution(sigma=1.5e6, mean_isotropic_chemical_shift=60.0)
    params = sf.make_LMFIT_params(spin_system_models=[model])

    calls = []
    original = sf._check_kernel_grid
    monkeypatch.setattr(
        sf, "_check_kernel_grid", lambda m, e: calls.append(1) or original(m, e)
    )

    for _ in range(3):
        sf._generate_distribution_spectrum(params, kernel, [model])
    assert len(calls) == 1

    # a new experiment on the same grid is checked again, once
    method.experiment = _experiment("27Al", 9.4, 256, 100e3, True, 1)
    for _ in range(2):
        sf._generate_distribution_spectrum(params, kernel, [model])
    assert len(calls) == 2

    # a mismatched experiment is still caught
    method.experiment = _experiment("27Al", 9.4, 256, 100e3, False, -1)
    with pytest.raises(ValueError, match="get_spectral_dimensions"):
        sf._generate_distribution_spectrum(params, kernel, [model])
