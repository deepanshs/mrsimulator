"""Anti-Phase test for an INADEQUATE method"""
import mrsimulator.signal_processor as sp
import numpy as np
from mrsimulator import Coupling
from mrsimulator import Method
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method import DelayEvent
from mrsimulator.method import MixingEvent
from mrsimulator.method import SpectralDimension
from mrsimulator.method import SpectralEvent


def setup_simulator_inadequate():
    j_coupling = 20
    site_A = Site(isotope="1H", isotropic_chemical_shift=-4.0)
    site_B = Site(isotope="1H", isotropic_chemical_shift=1.5)
    coupling_AB = Coupling(site_index=[0, 1], isotropic_j=j_coupling)
    sys = SpinSystem(sites=[site_A, site_B], couplings=[coupling_AB])

    delay = 1 / (4 * j_coupling)
    inadequate = Method(
        channels=["1H"],
        spectral_dimensions=[
            SpectralDimension(
                count=1024,
                spectral_width=8.0e3,
                reference_offset=0,
                label="2Q frequency",
                events=[
                    DelayEvent(
                        duration=2 * delay,
                        freq_contrib=["J"],
                        transition_queries=[{"ch1": {"P": [-1]}}],
                    ),
                    MixingEvent(query={"ch1": {"angle": np.pi / 2, "phase": 0.0}}),
                    SpectralEvent(
                        fraction=1.0, transition_queries=[{"ch1": {"P": [-1, -1]}}]
                    ),
                ],
            ),
            SpectralDimension(
                count=1024,
                spectral_width=4.0e3,
                reference_offset=0,
                label="1Q frequency",
                events=[
                    MixingEvent(
                        query={"ch1": {"angle": np.pi / 2, "phase": np.pi / 2}}
                    ),
                    SpectralEvent(
                        fraction=1.0, transition_queries=[{"ch1": {"P": [-1]}}]
                    ),
                ],
            ),
        ],
    )

    sim = Simulator(spin_systems=[sys], methods=[inadequate])
    sim.config.integration_density = 1
    sim.config.number_of_sidebands = 1

    return sim


def test_inadequate():
    sim = setup_simulator_inadequate()
    sim.run()

    # Apply processing to reduce intensity artifacts from binning
    operations = [
        sp.IFFT(dim_index=[0, 1]),
        sp.apodization.Gaussian(FWHM="100 Hz", dim_index=[0, 1]),
        sp.FFT(dim_index=[0, 1]),
    ]
    processor = sp.SignalProcessor(operations=operations)
    processed_data = processor.apply_operations(dataset=sim.methods[0].simulation).real
    processed_data /= processed_data.max()

    max_intensity = processed_data.max()
    min_intensity = processed_data.min()

    # Ensure peaks are anti phase
    assert np.allclose(max_intensity, -min_intensity, rtol=1e-4)
