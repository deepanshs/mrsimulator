import numpy as np
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method import Method
from mrsimulator.method import MixingEvent
from mrsimulator.method import SpectralDimension
from mrsimulator.method import SpectralEvent
from mrsimulator.method.query import MixingEnum
from mrsimulator.transition import TransitionPathway

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

method1 = Method(
    channels=["13C"],
    spectral_dimensions=[{"events": [{"transition_queries": [{"ch1": {"P": [-1]}}]}]}],
)

method2 = Method(
    channels=["13C"],
    spectral_dimensions=[{"events": [{"transition_queries": [{"ch1": {"P": [-2]}}]}]}],
)


def check_transition_set(got, expected):
    assert len(got) == len(expected), "Inconsistent transition pathway count"
    for item in got:
        assert item in expected, f"Transition pathways not found: {item}"


def test_00():
    system = SpinSystem(sites=[{"isotope": "13C"}])
    tr = method1.get_transition_pathways(system)

    expected = [TransitionPathway([{"final": [-0.5], "initial": [0.5]}])]
    check_transition_set(tr, expected)


def test_01():
    system = SpinSystem(sites=[{"isotope": "13C"}, {"isotope": "1H"}])
    tr = method1.get_transition_pathways(system)
    expected = [
        TransitionPathway([{"final": [-0.5, -0.5], "initial": [0.5, -0.5]}]),
        TransitionPathway([{"final": [-0.5, 0.5], "initial": [0.5, 0.5]}]),
    ]
    check_transition_set(tr, expected)


def test_02():
    system = SpinSystem(sites=[{"isotope": "13C"}, {"isotope": "13C"}])
    tr = method1.get_transition_pathways(system)

    expected = [
        TransitionPathway([{"final": [-0.5, -0.5], "initial": [0.5, -0.5]}]),
        TransitionPathway([{"final": [-0.5, 0.5], "initial": [0.5, 0.5]}]),
        TransitionPathway([{"final": [-0.5, -0.5], "initial": [-0.5, 0.5]}]),
        TransitionPathway([{"final": [0.5, -0.5], "initial": [0.5, 0.5]}]),
    ]
    check_transition_set(tr, expected)


def test_03():
    system = SpinSystem(sites=[{"isotope": "13C"}])
    tr = method2.get_transition_pathways(system)

    expected = []
    check_transition_set(tr, expected)


def test_04():
    system = SpinSystem(
        sites=[{"isotope": "13C"}, {"isotope": "13C"}, {"isotope": "14N"}]
    )
    tr = method1.get_transition_pathways(system)

    expected = [
        TransitionPathway([{"final": [-0.5, -0.5, -1], "initial": [0.5, -0.5, -1]}]),
        TransitionPathway([{"final": [-0.5, 0.5, -1], "initial": [0.5, 0.5, -1]}]),
        TransitionPathway([{"final": [-0.5, -0.5, -1], "initial": [-0.5, 0.5, -1]}]),
        TransitionPathway([{"final": [0.5, -0.5, -1], "initial": [0.5, 0.5, -1]}]),
        TransitionPathway([{"final": [-0.5, -0.5, 0], "initial": [0.5, -0.5, 0]}]),
        TransitionPathway([{"final": [-0.5, 0.5, 0], "initial": [0.5, 0.5, 0]}]),
        TransitionPathway([{"final": [-0.5, -0.5, 0], "initial": [-0.5, 0.5, 0]}]),
        TransitionPathway([{"final": [0.5, -0.5, 0], "initial": [0.5, 0.5, 0]}]),
        TransitionPathway([{"final": [-0.5, -0.5, 1], "initial": [0.5, -0.5, 1]}]),
        TransitionPathway([{"final": [-0.5, 0.5, 1], "initial": [0.5, 0.5, 1]}]),
        TransitionPathway([{"final": [-0.5, -0.5, 1], "initial": [-0.5, 0.5, 1]}]),
        TransitionPathway([{"final": [0.5, -0.5, 1], "initial": [0.5, 0.5, 1]}]),
    ]
    check_transition_set(tr, expected)


def test_hahn():
    system = SpinSystem(sites=[{"isotope": "13C"}, {"isotope": "13C"}])
    hahn = Method(
        channels=["13C"],
        spectral_dimensions=[
            {
                "events": [
                    {"fraction": 0.5, "transition_queries": [{"ch1": {"P": [1]}}]},
                    {"query": {"ch1": {"angle": np.pi, "phase": 0}}},
                    {"fraction": 0.5, "transition_queries": [{"ch1": {"P": [-1]}}]},
                ]
            },
        ],
    )
    tr = hahn.get_transition_pathways(system)

    weights = np.asarray([1, 1, 1, 1])
    transition_pathways = 0.5 * np.asarray(
        [
            [[[-1, -1], [-1, 1]], [[1, 1], [1, -1]]],
            [[[-1, -1], [1, -1]], [[1, 1], [-1, 1]]],
            [[[-1, 1], [1, 1]], [[1, -1], [-1, -1]]],
            [[[1, -1], [1, 1]], [[-1, 1], [-1, -1]]],
        ]
    )
    assert_transitions(transition_pathways, weights, tr)


def test_cosy():
    system = SpinSystem(sites=[{"isotope": "1H"}, {"isotope": "1H"}])
    cosy = Method(
        channels=["1H"],
        spectral_dimensions=[
            {
                "events": [
                    {"fraction": 1, "transition_queries": [{"ch1": {"P": [-1]}}]},
                    {"query": {"ch1": {"angle": np.pi / 2, "phase": 0}}},
                ],
            },
            {
                "events": [
                    {"fraction": 1, "transition_queries": [{"ch1": {"P": [-1]}}]}
                ],
            },
        ],
    )
    tr = cosy.get_transition_pathways(system)

    weights = np.asarray([1, -1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1, 1, 1, -1, 1]) * 0.25
    transition_pathways = 0.5 * np.asarray(
        [
            [[[-1, 1], [-1, -1]], [[-1, 1], [-1, -1]]],
            [[[-1, 1], [-1, -1]], [[1, -1], [-1, -1]]],
            [[[-1, 1], [-1, -1]], [[1, 1], [-1, 1]]],
            [[[-1, 1], [-1, -1]], [[1, 1], [1, -1]]],
            #
            [[[1, -1], [-1, -1]], [[-1, 1], [-1, -1]]],
            [[[1, -1], [-1, -1]], [[1, -1], [-1, -1]]],
            [[[1, -1], [-1, -1]], [[1, 1], [-1, 1]]],
            [[[1, -1], [-1, -1]], [[1, 1], [1, -1]]],
            #
            [[[1, 1], [-1, 1]], [[-1, 1], [-1, -1]]],
            [[[1, 1], [-1, 1]], [[1, -1], [-1, -1]]],
            [[[1, 1], [-1, 1]], [[1, 1], [-1, 1]]],
            [[[1, 1], [-1, 1]], [[1, 1], [1, -1]]],
            #
            [[[1, 1], [1, -1]], [[-1, 1], [-1, -1]]],
            [[[1, 1], [1, -1]], [[1, -1], [-1, -1]]],
            [[[1, 1], [1, -1]], [[1, 1], [-1, 1]]],
            [[[1, 1], [1, -1]], [[1, 1], [1, -1]]],
        ]
    )
    assert_transitions(transition_pathways, weights, tr)


def test_total_mixing():
    system = SpinSystem(sites=[{"isotope": "1H"}, {"isotope": "14N"}])
    total_mix = Method(
        channels=["1H"],
        spectral_dimensions=[
            SpectralDimension(
                events=[
                    SpectralEvent(transition_queries=[{"ch1": {"P": [-1]}}]),
                    MixingEvent(query=MixingEnum.TotalMixing),
                ]
            ),
            SpectralDimension(
                events=[SpectralEvent(transition_queries=[{"ch1": {"P": [-1]}}])]
            ),
        ],
    )
    transitions = total_mix.get_transition_pathways(system)
    tr_should_be = 0.5 * np.asarray(
        [
            [[[1, -2], [-1, -2]], [[1, -2], [-1, -2]]],
            [[[1, -2], [-1, -2]], [[1, 0], [-1, 0]]],
            [[[1, -2], [-1, -2]], [[1, 2], [-1, 2]]],
            [[[1, 0], [-1, 0]], [[1, -2], [-1, -2]]],
            [[[1, 0], [-1, 0]], [[1, 0], [-1, 0]]],
            [[[1, 0], [-1, 0]], [[1, 2], [-1, 2]]],
            [[[1, 2], [-1, 2]], [[1, -2], [-1, -2]]],
            [[[1, 2], [-1, 2]], [[1, 0], [-1, 0]]],
            [[[1, 2], [-1, 2]], [[1, 2], [-1, 2]]],
        ]
    )
    weights_should_be = np.ones(9)
    assert_transitions(tr_should_be, weights_should_be, transitions)


def test_no_mixing():
    system = SpinSystem(sites=[{"isotope": "1H"}, {"isotope": "14N"}])
    no_mix = Method(
        channels=["1H"],
        spectral_dimensions=[
            SpectralDimension(
                events=[
                    SpectralEvent(transition_queries=[{"ch1": {"P": [-1]}}]),
                    MixingEvent(query=MixingEnum.NoMixing),
                ]
            ),
            SpectralDimension(
                events=[SpectralEvent(transition_queries=[{"ch1": {"P": [-1]}}])]
            ),
        ],
    )
    transitions = no_mix.get_transition_pathways(system)
    tr_should_be = 0.5 * np.asarray(
        [
            [[[1, -2], [-1, -2]], [[1, -2], [-1, -2]]],
            [[[1, 0], [-1, 0]], [[1, 0], [-1, 0]]],
            [[[1, 2], [-1, 2]], [[1, 2], [-1, 2]]],
        ]
    )
    weights_should_be = np.ones(3)
    assert_transitions(tr_should_be, weights_should_be, transitions)


def assert_transitions(transition_pathways, weights, tr):
    expected = [
        TransitionPathway(
            pathway=[
                {"initial": list(states[0]), "final": list(states[1])}
                for states in transitions
            ],
            weight=w,
        )
        for transitions, w in zip(transition_pathways, weights)
    ]

    assert tr == expected


def test_empty_segments():
    sys = SpinSystem(sites=[Site(isotope="1H")])
    mth = Method(
        channels=["1H"],
        spectral_dimensions=[
            SpectralDimension(
                events=[
                    SpectralEvent(transition_queries=[{"ch1": {"P": [3]}}]),
                ]
            )
        ],
    )

    # Filter is empty since no p=3 transition for single proton, ensure no error raised
    assert mth._get_transition_pathways_np(sys).size == 0  # Empty list


def test_no_zero_transition_weights():
    """This method and spin system produces four transition pathways, two with weight
    1 and two with weight 0. Ensure 0 weight pathways are excluded
    """
    sys = SpinSystem(sites=[Site(isotope="2H")])
    mth = Method(
        channels=["2H"],
        spectral_dimensions=[
            SpectralDimension(
                events=[
                    SpectralEvent(transition_queries=[{"ch1": {"P": [1]}}]),
                    MixingEvent(query={"ch1": {"angle": np.pi}}),
                    SpectralEvent(transition_queries=[{"ch1": {"P": [-1]}}]),
                ]
            )
        ],
    )
    pathways = mth._get_transition_pathways_np(sys)
    weights = mth._get_transition_pathway_weights_np(pathways, sys)

    # Make sure zero pathways are present in the underlying calculation
    underlying_weights_should_be = np.asarray([0, 1, 1, 0])
    assert np.allclose(underlying_weights_should_be, weights)

    # Now, ensure zero weight pathways excluded with the combined pathway/weight func
    pathways, weights = mth._get_transition_pathway_and_weights_np(sys)

    transitions_should_be = np.asarray(
        [[[[-1], [0]], [[1], [0]]], [[[0], [1]], [[0], [-1]]]]
    )

    weights_should_be = np.ones(2)

    assert np.allclose(pathways, transitions_should_be)
    assert np.allclose(weights, weights_should_be)
