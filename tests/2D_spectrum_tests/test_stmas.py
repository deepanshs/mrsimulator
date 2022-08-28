import numpy as np
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method import Method
from mrsimulator.method import SpectralDimension
from mrsimulator.method.event import SpectralEvent
from mrsimulator.method.lib import ST1_VAS


def test_ST_VAS_isotropic_shift():
    Co59 = Site(
        isotope="59Co",
        isotropic_chemical_shift=-27.4,  # in ppm
        quadrupolar=dict(Cq=1.68e6, eta=0.2),  # Cq is in Hz
    )

    dim_kwargs = [
        dict(count=512, spectral_width=2e3, reference_offset=0),
        dict(count=512, spectral_width=1e3, reference_offset=-2e3),
    ]

    spin_systems = [SpinSystem(sites=[Co59])]
    method_2d = ST1_VAS(
        channels=["59Co"],
        magnetic_flux_density=7,  # in T
        rotor_angle=54.7359 * 3.14159 / 180,
        spectral_dimensions=[
            SpectralDimension(**dim_kwargs[0]),
            SpectralDimension(**dim_kwargs[1]),
        ],
    )
    method_1d_iso = Method(
        channels=["59Co"],
        magnetic_flux_density=7,  # in T
        rotor_angle=54.7359 * 3.14159 / 180,
        rotor_frequency=np.inf,
        spectral_dimensions=[
            SpectralDimension(
                **dim_kwargs[0],
                events=[
                    SpectralEvent(
                        fraction=1 / (1 + (84 / 135)),
                        freq_contrib=["Shielding1_0", "Quad2_0"],
                        transition_queries=[
                            {"ch1": {"P": [1], "D": [2]}},
                            {"ch1": {"P": [1], "D": [-2]}},
                        ],
                    ),
                    SpectralEvent(
                        fraction=(84 / 135) / (1 + (84 / 135)),
                        freq_contrib=["Shielding1_0", "Quad2_0"],
                        transition_queries=[{"ch1": {"P": [-1], "D": [0]}}],
                    ),
                ]
            ),
        ],
    )
    method_1d_ct = Method(
        channels=["59Co"],
        magnetic_flux_density=7,  # in T
        rotor_angle=54.7359 * 3.14159 / 180,
        rotor_frequency=np.inf,
        spectral_dimensions=[
            SpectralDimension(
                **dim_kwargs[1],
                events=[
                    SpectralEvent(transition_queries=[{"ch1": {"P": [-1], "D": [0]}}])
                ]
            ),
        ],
    )
    methods = [method_2d, method_1d_iso, method_1d_ct]

    sim = Simulator(spin_systems=spin_systems, methods=methods)
    sim.run()

    dataset_proj_iso = sim.methods[0].simulation.sum(axis=0)
    dataset_proj_ct = sim.methods[0].simulation.sum(axis=1)
    dataset_iso = sim.methods[1].simulation
    dataset_ct = sim.methods[2].simulation

    # check isotropic freq is the same
    m1 = np.argmax(dataset_proj_iso.y[0].components[0])
    m2 = np.argmax(dataset_iso.y[0].components[0])
    assert m1 == m2

    # check the quad CT spectrum is same.
    ct_proj = dataset_proj_ct.y[0].components[0]
    ct_spec = dataset_ct.y[0].components[0]
    np.testing.assert_almost_equal(
        ct_proj / ct_proj.sum(), ct_spec / ct_spec.sum(), decimal=8
    )
