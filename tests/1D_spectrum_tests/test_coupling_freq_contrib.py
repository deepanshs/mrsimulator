import numpy as np
from mrsimulator import Coupling
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method import Method
from mrsimulator.method import MixingEventA
from mrsimulator.method import SpectralDimension
from mrsimulator.method import SpectralEvent
from mrsimulator.spin_system.tensors import SymmetricTensor

positive_sq_tq = [{"ch1": {"P": [1]}}]
negative_sq_tq = [{"ch1": {"P": [-1]}}]


def coupled_spin_system(j_coup=0, dipole=0):
    S1 = Site(
        isotope="1H",
        isotropic_chemical_shift=10,  # in ppm
        shielding_symmetric=SymmetricTensor(zeta=-80, eta=0.25),  # zeta in ppm
    )
    S2 = Site(isotope="1H", isotropic_chemical_shift=-10)
    S12 = Coupling(
        site_index=[0, 1],
        isotropic_j=j_coup,
        dipolar=SymmetricTensor(D=dipole, eta=0),
    )
    return SpinSystem(sites=[S1, S2], couplings=[S12])


def hahn_method():
    return Method(
        channels=["1H"],
        magnetic_flux_density=9.4,  # in T
        rotor_angle=0,  # in rads
        spectral_dimensions=[
            SpectralDimension(
                count=512,
                spectral_width=2e4,  # in Hz
                events=[
                    SpectralEvent(fraction=0.5, transition_queries=positive_sq_tq),
                    MixingEventA(ch1={"angle": np.pi, "phase": 0}),
                    SpectralEvent(fraction=0.5, transition_queries=negative_sq_tq),
                ],
            )
        ],
    )


def simulator(method, j_coup, dipole):
    system = coupled_spin_system(j_coup, dipole)
    sim = Simulator(spin_systems=[system], methods=[method])
    sim.run()
    return sim.methods[0].simulation.real.y[0].components[0]


def contrib_method(contrib):
    return Method(
        channels=["1H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=512,
                spectral_width=2e4,
                events=[
                    SpectralEvent(
                        freq_contrib=contrib, transition_queries=negative_sq_tq
                    )
                ],
            )
        ],
    )


def test_01():
    hahn = hahn_method()
    hahn_sim = simulator(hahn, j_coup=500, dipole=1000)

    test_method = contrib_method(["J1_0", "D1_2"])
    test_sim = simulator(test_method, j_coup=500, dipole=1000)

    np.testing.assert_almost_equal(hahn_sim, test_sim, decimal=8)


def test_02():
    hahn = hahn_method()
    hahn_sim = simulator(hahn, j_coup=0, dipole=1000)

    test_method = contrib_method(["D1_2"])
    test_sim = simulator(test_method, j_coup=500, dipole=1000)

    np.testing.assert_almost_equal(hahn_sim, test_sim, decimal=8)


def test_03():
    hahn = hahn_method()
    hahn_sim = simulator(hahn, j_coup=5000, dipole=0)

    test_method = contrib_method(["J1_0"])
    test_sim = simulator(test_method, j_coup=5000, dipole=2000)

    np.testing.assert_almost_equal(hahn_sim, test_sim, decimal=8)
