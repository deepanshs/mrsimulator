# -*- coding: utf-8 -*-
import numpy as np
import pytest
from mrsimulator.method import Method
from mrsimulator.method.plot import _add_tip_angle_and_phase
from mrsimulator.method.plot import _check_columns
from mrsimulator.method.plot import _format_df
from mrsimulator.method.plot import _make_normal_and_offset_x_data
from mrsimulator.method.plot import _make_x_data
from mrsimulator.method.plot import _offset_x_data

__author__ = "Matthew D. Giammar"
__email__ = "giammar.7@buckeyemail.osu.edu"


def test_check_columns():
    df1 = method1_df()
    df2 = method2_df()

    # Test required columns present in DataFrame
    bad_df1 = df1.drop("fraction", axis=1)
    bad_df2 = df2.drop(["p", "d"], axis=1)
    error = r".*Some required columns were not present in the DataFrame.*"
    with pytest.raises(ValueError, match=error):
        _check_columns(bad_df1)
    with pytest.raises(ValueError, match=error):
        _check_columns(bad_df2)

    # Test correct columns dropped from DataFrame
    assert _check_columns(df1) == []
    assert _check_columns(df2) == [
        "magnetic_flux_density",
        "rotor_frequency",
        "rotor_angle",
    ]


def test_add_tip_angle_and_phase():
    df1 = method1_df()
    df2 = method2_df()

    # Modify the dataframes
    _add_tip_angle_and_phase(df1)
    _add_tip_angle_and_phase(df2)

    # Check 'tip_angle' and 'phase' are columns now in DataFrame
    assert "tip_angle" in df1.columns
    assert "phase" in df1.columns
    assert "tip_angle" in df2.columns
    assert "phase" in df2.columns

    # Check correct calculations for tip angle and phase
    ta_should_be = np.array([np.nan, np.nan, 90.0, np.nan])
    p_should_be = np.array([np.nan, np.nan, 0.0, np.nan])

    assert np.allclose(np.asarray(df1["tip_angle"]), ta_should_be, equal_nan=True)
    assert np.allclose(np.asarray(df1["phase"]), p_should_be, equal_nan=True)

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

    assert np.allclose(np.asarray(df2["tip_angle"]), ta_should_be, equal_nan=True)
    assert np.allclose(np.asarray(df2["phase"]), p_should_be, equal_nan=True)


def test_format_df():
    df1 = method1_df()
    df2 = method2_df()
    params1 = _format_df(df1)
    params2 = _format_df(df2)

    # Check 'tip_angle' and 'phase' are columns now in DataFrame
    assert "tip_angle" in df1.columns
    assert "phase" in df1.columns
    assert "tip_angle" in df2.columns
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

    # Check expected x_data retunred
    # NOTE: Should arrays be hardcoded? Or should be calculated in simmilar way
    x1_should_be = [0, 0.8, 0.8, 1.3, 1.3, 1.7]
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
    _add_tip_angle_and_phase(df1)
    _add_tip_angle_and_phase(df2)

    # Calculate needed x_data
    x1 = _make_x_data(df1)
    x2 = _make_x_data(df2)
    off_x1 = _offset_x_data(df1, x1)
    off_x2 = _offset_x_data(df2, x2)

    # Check return type
    assert isinstance(off_x1, np.ndarray)
    assert isinstance(off_x2, np.ndarray)

    # Check expected offset_x returned
    # NOTE: Should arrays be hardcoded? Or should be calculated in simmilar way
    off_x1_should_be = [0.0, 0.0, 0.8, 0.8, 1.26875, 1.33125, 1.7]
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

    assert np.allclose(off_x1, off_x1_should_be)
    assert np.allclose(off_x2, off_x2_should_be)


def test_make_normal_and_offset_x_data():
    # Setup DataFrames
    df1 = method1_df()
    df2 = method2_df()
    _add_tip_angle_and_phase(df1)
    _add_tip_angle_and_phase(df2)

    # Test when x_data is lenght zero
    no_events_df1 = df1.drop(df1.index, axis=0)
    no_events_df2 = df2.drop(df2.index[1:], axis=0)

    error = (
        r".*The DataFrame does not contain any SpectralEvents or "
        r"ConstandDurationEvents. At least one must be present to construct a plot.*"
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
    assert isinstance(off_x1, np.ndarray)
    assert isinstance(off_x2, np.ndarray)


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
                    {"mixing_query": {"ch1": {"tip_angle": np.pi / 2, "phase": 0}}},
                    {"fraction": 0.5},
                ],
            },
        ],
    )

    return method1.summary()


def method2_df():
    method2 = Method(
        channels=["1H"],
        spectral_dimensions=[
            {
                "label": "Mixing, Spectral, Mixing, Mixing, ConstantDuration ",
                "events": [
                    {"mixing_query": {"ch1": {"tip_angle": np.pi / 2, "phase": 0.3}}},
                    {"fraction": 1},
                    {"mixing_query": {"ch1": {"tip_angle": np.pi, "phase": 0.15}}},
                    {"mixing_query": {"ch1": {"tip_angle": 3 * np.pi / 2, "phase": 0}}},
                    {"duration": 1.5},
                ],
            },
            {
                "label": "Spectral, Spectral, CondatntDuration, ConstantDuration",
                "events": [
                    {"fraction": 0.2},
                    {"fraction": 0.8},
                    {"duration": 1.3},
                    {"duration": 0.6},
                ],
            },
            {
                "label": "Mixing, Mixing, Mixing, Spectral, ConstantDuration, Spectal",
                "events": [
                    {"mixing_query": {"ch1": {"tip_angle": 1.2, "phase": np.pi}}},
                    {"mixing_query": {"ch1": {"tip_angle": 0.75, "phase": np.pi / 4}}},
                    {"mixing_query": {"ch1": {"tip_angle": 0.33, "phase": np.pi / 3}}},
                    {"fraction": 0.5},
                    {"duration": 1.5},
                    {"fraction": 0.5},
                ],
            },
        ],
    )

    return method2.summary(drop_constant_columns=False)
