import mrsimulator.signal_processor as sp
import numpy as np
from mrsimulator import Coupling
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method import DelayEvent
from mrsimulator.method import Method
from mrsimulator.method import MixingEventA
from mrsimulator.method import SpectralDimension
from mrsimulator.method import SpectralEvent


tq_plus_one = {"ch1": {"P": [+1]}}
tq_minus_one = {"ch1": {"P": [-1]}}
tq_2Q = {"ch1": {"P": [-1, -1]}}


def build_2Q_method_full(delay):
    """Build a 2Q method object with a given delay in seconds"""
    return Method(
        channels=["1H"],
        spectral_dimensions=[
            SpectralDimension(
                count=8192,
                spectral_width=2.5e3,
                reference_offset=0,
                events=[
                    DelayEvent(duration=delay, transition_queries=[tq_plus_one]),
                    MixingEventA(ch1={"angle": np.pi, "phase": 0}),
                    DelayEvent(duration=delay, transition_queries=[tq_minus_one]),
                    MixingEventA(ch1={"angle": np.pi / 2, "phase": 0}),
                    SpectralEvent(fraction=1.0, transition_queries=[tq_2Q]),
                ],
            )
        ],
    )


def build_2Q_method_simplified(delay):
    """Build a 2Q method object where the refocusing has been simplified with a given
    delay in seconds.
    """
    return Method(
        channels=["1H"],
        spectral_dimensions=[
            SpectralDimension(
                count=8192,
                spectral_width=2.5e3,
                reference_offset=0,
                events=[
                    DelayEvent(
                        duration=2 * delay,
                        transition_queries=[tq_minus_one],
                        freq_contrib=["J"],
                    ),
                    MixingEventA(ch1={"angle": np.pi / 2, "phase": 0}),
                    SpectralEvent(fraction=1.0, transition_queries=[tq_2Q]),
                ],
            )
        ],
    )


def build_J_coupled_system(j_coupling):
    """Create a coupled spin system with a given J coupling strength"""
    sites = np.asarray(
        [
            Site(isotope="1H", isotropic_chemical_shift=-2.0),  # in ppm
            Site(isotope="1H", isotropic_chemical_shift=1.5),  # in ppm
        ]
    )
    couplings = np.asarray([Coupling(site_index=[0, 1], isotropic_j=j_coupling)])
    return SpinSystem(sites=sites, couplings=couplings)


def setup_and_process_simulation(sys, methods):
    # Create and run simulator with single orientation (liquid simulation)
    sim = Simulator(spin_systems=[sys], methods=methods)
    sim.config.integration_density = 1
    sim.config.number_of_sidebands = 1
    sim.run()

    operations = [
        sp.IFFT(dim_index=0),
        sp.apodization.Gaussian(FWHM="4 Hz", dim_index=0),
        sp.FFT(dim_index=0),
    ]
    processor = sp.SignalProcessor(operations=operations)
    return [processor.apply_operations(dataset=mth.simulation) for mth in sim.methods]


def test_2Q_buildup():
    j_range = np.asarray([1.0, 2.5, 5.0, 10.0, 20.0, 50.0, 100.0])

    # Iterate over different J coupling strengths and assert maximum happens at tau=1/4J
    for j_coupling in j_range:
        sys = build_J_coupled_system(j_coupling)

        # Tau ranges from 0.0 to 1/2J with 1/4J at center
        tau_range = np.linspace(start=0.0, stop=2.0, num=31) / (4 * j_coupling)
        methods_full = [build_2Q_method_full(delay) for delay in tau_range]

        # Get processed spectra and find maximum intensities
        processed_data = setup_and_process_simulation(sys, methods_full)
        max_intensities = [data.max().value for data in processed_data]

        # Check that maximum intensity for methods happens at delay = 1/4J
        assert (
            np.argmax(max_intensities) == np.where(tau_range == 1 / (4 * j_coupling))[0]
        )

        # Same test, but for simplified method using freq_contrib
        methods_simplified = [build_2Q_method_simplified(delay) for delay in tau_range]

        # Get processed spectra and find maximum intensities
        processed_data = setup_and_process_simulation(sys, methods_simplified)
        max_intensities = [data.max().value for data in processed_data]

        # Check that maximum intensity for methods happens at delay = 1/4J
        assert (
            np.argmax(max_intensities) == np.where(tau_range == 1 / (4 * j_coupling))[0]
        )
