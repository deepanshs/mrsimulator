import matplotlib.projections as proj
import numpy as np
import pandas as pd
import pytest
from matplotlib import pyplot as plt
from mrsimulator.method import Method
from mrsimulator.method.plot import _add_angle_and_phase
from mrsimulator.method.plot import _check_columns
from mrsimulator.method.plot import _format_df
from mrsimulator.method.plot import _make_normal_and_offset_x_data
from mrsimulator.method.plot import _make_x_data
from mrsimulator.method.plot import _offset_x_data
from mrsimulator.method.plot import CustomAxes
from mrsimulator.method.plot import MultiLineAxes
from mrsimulator.method.plot import SequenceDiagram

# from mrsimulator.method.plot import SpectralDimensionLabels


__author__ = "Matthew D. Giammar"
__email__ = "giammar.7@buckeyemail.osu.edu"


def test_check_columns():
    df1 = method1_df()
    df2 = method2_df()

    # Test required columns present in DataFrame
    bad_df1 = df1.drop("fraction", axis=1)
    bad_df2 = df2.drop(["p"], axis=1)
    error = r".*Some required columns were not present in the DataFrame.*"
    with pytest.raises(ValueError, match=error):
        _check_columns(bad_df1)
    with pytest.raises(ValueError, match=error):
        _check_columns(bad_df2)

    # Test correct columns dropped from DataFrame
    assert _check_columns(df1) == ["d"]
    assert _check_columns(df2) == [
        "magnetic_flux_density",
        "rotor_frequency",
        "rotor_angle",
        "d",
    ]


def test_add_angle_and_phase():
    df1 = method1_df()
    df2 = method2_df()

    # Modify the dataframes
    _add_angle_and_phase(df1)
    _add_angle_and_phase(df2)

    # Check 'angle' and 'phase' are columns now in DataFrame
    assert "angle" in df1.columns
    assert "phase" in df1.columns
    assert "angle" in df2.columns
    assert "phase" in df2.columns

    # Check correct calculations for tip angle and phase
    ta_should_be = np.array([np.nan, np.nan, np.nan, 90.0, np.nan])
    p_should_be = np.array([np.nan, np.nan, np.nan, 0.0, np.nan])

    df1_ch1_tip = [x["ch1"] if x is not None else np.nan for x in df1["angle"]]
    df1_ch1_phase = [x["ch1"] if x is not None else np.nan for x in df1["phase"]]

    assert np.allclose(np.asarray(df1_ch1_tip), ta_should_be, equal_nan=True)
    assert np.allclose(np.asarray(df1_ch1_phase), p_should_be, equal_nan=True)

    ta_should_be = np.array(
        [
            90.0,
            np.nan,
            180.0,
            270.0,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            68.75493541569878,
            42.97183463481174,
            18.90760723931717,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    p_should_be = np.array(
        [
            17.188733853924695,
            np.nan,
            8.594366926962348,
            0.0,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            180.0,
            45.0,
            60.0,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    ta_should_be_ch2 = np.array(
        [
            10,
            np.nan,
            20,
            30,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            40,
            50,
            60,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    p_should_be_ch2 = ta_should_be_ch2

    df2_ch1_tip = [x["ch1"] if x is not None else np.nan for x in df2["angle"]]
    df2_ch2_tip = [x["ch2"] if x is not None else np.nan for x in df2["angle"]]
    df2_ch1_phase = [x["ch1"] if x is not None else np.nan for x in df2["phase"]]
    df2_ch2_phase = [x["ch2"] if x is not None else np.nan for x in df2["phase"]]

    assert np.allclose(np.asarray(df2_ch1_tip), ta_should_be, equal_nan=True)
    assert np.allclose(np.asarray(df2_ch2_tip), ta_should_be_ch2, equal_nan=True)

    assert np.allclose(np.asarray(df2_ch1_phase), p_should_be, equal_nan=True)
    assert np.allclose(np.asarray(df2_ch2_phase), p_should_be_ch2, equal_nan=True)


def test_format_df():
    df1 = method1_df()
    df2 = method2_df()
    params1 = _format_df(df1)
    params2 = _format_df(df2)

    # 'd' should not be in params1 since 'd' is never defined
    assert "d" not in params1

    # Check 'angle' and 'phase' are columns now in DataFrame
    assert "angle" in df1.columns
    assert "phase" in df1.columns
    assert "angle" in df2.columns
    assert "phase" in df2.columns

    # Check expected params returned
    assert isinstance(params1, list)
    assert isinstance(params2, list)


def test_make_x_data():
    df1 = method1_df()
    df2 = method2_df()
    x1 = _make_x_data(df1)
    x2 = _make_x_data(df2)

    # Check return type
    assert isinstance(x1, list)
    assert isinstance(x2, list)

    # Check expected x_data returned
    x1_should_be = [0, 0.8, 0.8, 1.3, 1.3, 1.7, 1.7, 2.1]
    x2_should_be = [
        0,
        0.8,
        0.8,
        1.3,
        1.3,
        1.46,
        1.46,
        2.1,
        2.1,
        2.6,
        2.6,
        3.1,
        3.1,
        3.5,
        3.5,
        4.0,
        4.0,
        4.4,
    ]

    assert np.allclose(x1, x1_should_be)
    assert np.allclose(x2, x2_should_be)


def test_offset_x_data():
    # Setup DataFrames
    df1 = method1_df()
    df2 = method2_df()
    df3 = method3_df()
    _add_angle_and_phase(df1)
    _add_angle_and_phase(df2)
    _add_angle_and_phase(df3)

    # Calculate needed x_data
    x1 = _make_x_data(df1)
    x2 = _make_x_data(df2)
    x3 = _make_x_data(df3)
    off_x1 = _offset_x_data(df1, x1)
    off_x2 = _offset_x_data(df2, x2)
    off_x3 = _offset_x_data(df3, x3)

    # Check return type
    # assert isinstance(off_x1, dict)
    # assert isinstance(off_x2, dict)
    # assert isinstance(off_x3, dict)
    assert list(off_x1.keys()) == ["ch1"]
    assert list(off_x2.keys()) == ["ch1", "ch2"]
    assert list(off_x3.keys()) == ["ch1"]

    for item in off_x1.values():
        assert isinstance(item, np.ndarray)
    for item in off_x2.values():
        assert isinstance(item, np.ndarray)
    for item in off_x3.values():
        assert isinstance(item, np.ndarray)

    # Check expected offset_x returned
    off_x1_should_be = [0.0, 0.0, 0.8, 0.8, 1.3, 1.3, 1.66875, 1.73125, 2.1]
    off_x2_should_be = [
        0.0,
        0.0625,
        0.64375,
        0.95625,
        1.3,
        1.3,
        1.46,
        1.46,
        2.1,
        2.1,
        2.6,
        2.6,
        3.0546408412188097,
        3.1453591587811904,
        3.5,
        3.5,
        4.0,
        4.0,
        4.4,
    ]
    off_x3_should_be = [0, 0, 0.08, 0.08, 0.8]

    assert np.allclose(off_x1["ch1"], off_x1_should_be)
    assert np.allclose(off_x2["ch1"], off_x2_should_be)
    assert np.allclose(off_x3["ch1"], off_x3_should_be)


def test_make_normal_and_offset_x_data():
    # Setup DataFrames
    df1 = method1_df()
    df2 = method2_df()
    _add_angle_and_phase(df1)
    _add_angle_and_phase(df2)

    # Test when x_data is length zero
    no_events_df1 = df1.drop(df1.index, axis=0)
    no_events_df2 = df2.drop(df2.index[1:], axis=0)

    error = (
        r".*The DataFrame does not contain any SpectralEvents or "
        r"ConstantDurationEvents. At least one must be present to construct a plot.*"
    )
    with pytest.raises(ValueError, match=error):
        _make_normal_and_offset_x_data(no_events_df1)
    with pytest.raises(ValueError, match=error):
        _make_normal_and_offset_x_data(no_events_df2)

    # Check return types
    x1, off_x1 = _make_normal_and_offset_x_data(df1)
    x2, off_x2 = _make_normal_and_offset_x_data(df2)

    assert isinstance(x1, list)
    assert isinstance(x2, list)
    assert isinstance(off_x1, dict)
    assert isinstance(off_x2, dict)


def test_SequenceDiagram():
    df3 = method3_df()
    _add_angle_and_phase(df3)
    x3 = [0.0, 0.0, 0.08, 0.08, 0.8]
    x_data_should_be = [0, 0, 0.08, 0.08, 0.675, 0.8]

    # Setup matplotlib objects
    fig = plt.figure()
    proj.register_projection(SequenceDiagram)
    axes_obj = fig.add_subplot(projection="sequence_axes")
    axes_obj.plot_diagram(df3, x3, "1H", "ch1")

    assert np.allclose(axes_obj.x_data, x_data_should_be)


# def test_SpectralDimensionLabels():
#     df2 = method2_df()
#     _add_angle_and_phase(df2)
#     x2 = _make_x_data(df2)
#     off_x2 = _offset_x_data(df2, x2)

#     # Setup matplotlib objects
#     fig = plt.figure()
#     proj.register_projection(SpectralDimensionLabels)
#     axes_obj = fig.add_subplot(projection="label_axes")
#     axes_obj.plot_labels(df2, off_x2)


def test_MultilineAxes():
    x_data = [0, 0, 1, 1, 2, 2, 3, 3]
    y_data = pd.Series([[1, 1, 1], [1, 0, 0], [1, 0, -1], [1, 1, 0]])

    # Setup matplotlib objects
    fig = plt.figure()
    proj.register_projection(MultiLineAxes)
    axes_obj = fig.add_subplot(projection="multi_line_axes")

    # make_plot
    # Test y_data dimensions check
    bad_y_data = pd.Series([0, 1, 2, 3])  # not 2 dimensional
    error = r".*Symmetry pathway data is misshapen. Data must be 2d.*"
    with pytest.raises(ValueError, match=error):
        axes_obj.make_plot(x_data, bad_y_data, "", [], {}, {})

    # _offset_overlaps
    # Test only one symmetry pathway
    y_data = np.zeros((1, 10))
    assert np.array_equal(axes_obj._offset_overlaps(y_data), y_data)

    # Test overlapping symmetry pathways
    # Each column represents a symmetry pathway
    y_data = np.array(
        [
            np.array([1, 1, 1, 1], dtype=float),
            np.array([1, 0, 0, 1], dtype=float),
            np.array([1, 0, -1, 0], dtype=float),
        ]
    )
    offset_should_be = [
        [0.91, 0.97, 1.0, 0.94],
        [1.0, -0.06, 0.0, 1.06],
        [1.09, 0.06, -1.0, -0.03],
    ]
    assert np.allclose(axes_obj._offset_overlaps(y_data), offset_should_be)


def test_CustomAxes():
    # Setup matplotlib objects
    fig = plt.figure()
    proj.register_projection(CustomAxes)
    axes_obj = fig.add_subplot(projection="custom_axes")

    # _add_rect_with_labels
    # Test error thrown when no color supplied
    error = r".*No color in `rect_kwargs`. A color must be specified.*"
    with pytest.raises(ValueError, match=error):
        axes_obj._add_rect_with_label(0, 1, "", {"height": 1})

    # _add_blank_space_labels
    axes_obj.make_plot(
        [1, 2, 2, 3, 3, 4], [0, 1e12, np.nan], "rotor_frequency", [False] * 3, {}, {}
    )


def method1_df():
    method1 = Method(
        channels=["1H"],
        spectral_dimensions=[
            {
                "label": "Spectral and Constant Duration Event",
                "events": [
                    {"fraction": 1},
                    {"duration": 1.5},
                ],
            },
            {
                "label": "Mixing and Spectral Event",
                "events": [
                    {"fraction": 0.5},
                    {"query": {"ch1": {"angle": np.pi / 2, "phase": 0}}},
                    {"fraction": 0.5},
                ],
            },
        ],
    )

    return method1.summary()


def method2_df():
    factor = np.pi / 180.0
    method2 = Method(
        channels=["1H"],
        spectral_dimensions=[
            {
                "label": "Mixing, Spectral, Mixing, Mixing, ConstantDuration ",
                "events": [
                    {
                        "query": {
                            "ch1": {"angle": np.pi / 2, "phase": 0.3},
                            "ch2": {"angle": 10 * factor, "phase": 10 * factor},
                        }
                    },
                    {
                        "fraction": 1,
                        "transition_query": [{"ch1": {"P": [1], "D": [0]}}],
                    },
                    {
                        "query": {
                            "ch1": {"angle": np.pi, "phase": 0.15},
                            "ch2": {"angle": 20 * factor, "phase": 20 * factor},
                        }
                    },
                    {
                        "query": {
                            "ch1": {"angle": 3 * np.pi / 2, "phase": 0},
                            "ch2": {"angle": 30 * factor, "phase": 30 * factor},
                        }
                    },
                    {
                        "duration": 1.5,
                        "transition_query": [{"ch1": {"P": [-1], "D": [2]}}],
                    },
                ],
            },
            {
                "label": "Spectral, Spectral, ConstantDuration, ConstantDuration",
                "events": [
                    {"fraction": 0.2},
                    {"fraction": 0.8},
                    {"duration": 1.3},
                    {"duration": 0.6},
                ],
            },
            {
                "label": "Mixing, Mixing, Mixing, Spectral, ConstantDuration, Spectral",
                "events": [
                    {
                        "query": {
                            "ch1": {"angle": 1.2, "phase": np.pi},
                            "ch2": {"angle": 40 * factor, "phase": 40 * factor},
                        }
                    },
                    {
                        "query": {
                            "ch1": {"angle": 0.75, "phase": np.pi / 4},
                            "ch2": {"angle": 50 * factor, "phase": 50 * factor},
                        }
                    },
                    {
                        "query": {
                            "ch1": {"angle": 0.33, "phase": np.pi / 3},
                            "ch2": {"angle": 60 * factor, "phase": 60 * factor},
                        }
                    },
                    {"fraction": 0.5},
                    {"duration": 1.5},
                    {"fraction": 0.5},
                ],
            },
        ],
    )

    return method2.summary(drop_constant_columns=False)


# something with calculating total tip angle with method 3
def method3_df():
    method3 = Method(
        channels=["1H"],
        spectral_dimensions=[
            {
                "label": "Spectral Spectral Mixing",
                "events": [
                    {
                        "fraction": 0.1,
                        "transition_query": [
                            {"ch1": {"P": [1, 1]}},
                            {"ch1": {"P": [2]}},
                        ],
                    },
                    {
                        "fraction": 0.9,
                        "transition_query": [
                            {"ch1": {"P": [-1, -1]}},
                            {"ch1": {"P": [-2]}},
                        ],
                    },
                    {"query": {"ch1": {"angle": np.pi, "phase": 0}}},
                ],
            }
        ],
    )

    return method3.summary()
