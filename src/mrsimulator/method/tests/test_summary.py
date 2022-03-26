# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from matplotlib.pyplot import Figure
from mrsimulator.method import Method

__author__ = "Matthew D. Giammar"
__email__ = "giammar.7@buckeyemail.osu.edu"


ME = "MixingEvent"
CDE = "ConstantDurationEvent"
SE = "SpectralEvent"
REQUIRED = [
    "type",
    "label",
    "duration",
    "fraction",
    "mixing_query",
    "spec_dim_index",
    "spec_dim_label",
    "p",
    "d",
]
ALL_PARAMS = [
    "type",
    "spec_dim_index",
    "spec_dim_label",
    "label",
    "duration",
    "fraction",
    "mixing_query",
    "magnetic_flux_density",
    "rotor_frequency",
    "rotor_angle",
    "freq_contrib",
    "p",
    "d",
]


def basic_summary_tests(the_method):
    def check_col_equal(col, arr):
        check = []
        for col_item, arr_item in zip(col, arr):
            if np.isnan(col_item):
                check.append(np.isnan(arr_item))
            else:
                check.append(np.isclose(col_item, arr_item))
        # Returns True if all items are equal, False otherwise
        return np.all(check)

    def check_col_equal_2d(col, arr):
        check = []
        for col_list, arr_list in zip(col, arr):
            for col_item, arr_item in zip(col_list, arr_list):
                if np.isnan(col_item):
                    check.append(np.isnan(arr_item))
                else:
                    check.append(np.isclose(col_item, arr_item))
        # Returns True if all items are equal, False otherwise
        return np.all(check)

    df = the_method.summary(drop_constant_columns=False)

    # Check correct return type
    assert isinstance(df, pd.DataFrame)

    # Check for correct number of event and properties (df.shape)
    assert df.shape[0] == 6
    assert df.shape[1] == len(ALL_PARAMS)

    # Check for correct columns
    assert set(df.columns) == set(ALL_PARAMS)

    # Check type
    assert np.all(df["type"] == [ME, CDE, SE, ME, CDE, SE])

    # Check label
    assert np.all(df["label"] == ["Mix0", "Dur0", "Spec0", "Mix1", "Dur1", "Spec1"])

    # NOTE: Need to run loops so NANs can be checked since `nan == nan` evals False
    # Check duration
    temp = [np.nan, 1.0, np.nan, np.nan, 2.0, np.nan]
    assert check_col_equal(df["duration"], temp)

    # Check fraction
    temp = [np.nan, np.nan, 1, np.nan, np.nan, 1]
    assert check_col_equal(df["fraction"], temp)

    # Check p
    temp = [[np.nan], [0.0], [1.0], [np.nan], [3.0], [4.0]]
    assert check_col_equal_2d(df["p"], temp)

    # Check d
    temp = [[np.nan], [0.0], [-2.0], [np.nan], [-4.0], [-6.0]]
    assert check_col_equal_2d(df["d"], temp)

    # Check magnetic_flux_density
    temp = [np.nan, 1.0, 2.0, np.nan, 3.0, 4.0]
    assert check_col_equal(df["magnetic_flux_density"], temp)

    # Check rotor_frequency
    temp = [np.nan, 0, 20.0, np.nan, 0, 0]
    assert check_col_equal(df["rotor_frequency"], temp)

    # Check rotor_angle (degrees)
    temp = [
        np.nan,
        5.729577951308232,
        11.459155902616464,
        np.nan,
        17.188733853924695,
        22.91831180523293,
    ]
    assert check_col_equal(df["rotor_angle"], temp)

    assert "mixing_query" in df.columns
    assert "freq_contrib" in df.columns


def args_summary_tests(the_method):
    # Pass drop_constant_columns=False
    df = the_method.summary(drop_constant_columns=False)

    assert df.shape[0] == 6
    assert df.shape[1] == len(ALL_PARAMS)
    assert set(df.columns) == set(ALL_PARAMS)

    # Pass drop_constant_columns=True (default)
    df = the_method.summary(drop_constant_columns=True)
    props_should_be = ["rotor_frequency", "freq_contrib"]

    assert df.shape[0] == 6
    assert df.shape[1] == len(REQUIRED + props_should_be)
    assert set(df.columns) == set(REQUIRED + props_should_be)

    # Make sure columns always present are not dropped even if constant
    event_ = [{"label": "all labels the same", "fraction": 0.2}]
    all_SpectralEvents_spec_dims = [{"events": event_ * 5}]

    the_method = Method(
        name="all-spectral-method-1",
        channels=["1H", "13C"],
        spectral_dimensions=all_SpectralEvents_spec_dims,
    )

    df = the_method.summary(drop_constant_columns=True)

    assert df.shape[0] == 5
    assert df.shape[1] == len(REQUIRED + ["freq_contrib"])
    assert set(df.columns) == set(REQUIRED + ["freq_contrib"])


def test_summary():
    all_defined_no_constant_spec_dims = [
        {
            "events": [
                {
                    "label": "Mix0",
                    "mixing_query": {
                        "ch1": {"tip_angle": np.pi / 4, "phase": np.pi / 2}
                    },
                },
                {
                    "label": "Dur0",
                    "duration": 1,
                    "magnetic_flux_density": 1,
                    "rotor_frequency": 0,  # in kHz
                    "rotor_angle": 0.1,
                    "transition_query": [{"ch1": {"P": [0], "D": [0]}}],
                },
                {
                    "label": "Spec0",
                    "fraction": 1,
                    "magnetic_flux_density": 2,
                    "rotor_frequency": 20000,
                    "rotor_angle": 0.2,
                    "transition_query": [{"ch1": {"P": [1], "D": [-2]}}],
                },
            ]
        },
        {
            "events": [
                {
                    "label": "Mix1",
                    "mixing_query": {"ch1": {"tip_angle": np.pi / 2, "phase": np.pi}},
                },
                {
                    "label": "Dur1",
                    "duration": 2,
                    "magnetic_flux_density": 3,
                    "rotor_frequency": 0,
                    "rotor_angle": 0.3,
                    "transition_query": [{"ch1": {"P": [3], "D": [-4]}}],
                },
                {
                    "label": "Spec1",
                    "fraction": 1,
                    "magnetic_flux_density": 4,
                    "rotor_frequency": 0,
                    "rotor_angle": 0.4,
                    "transition_query": [{"ch1": {"P": [4], "D": [-6]}}],
                },
            ]
        },
    ]

    const_mfd_rotor_ang_spec_dims = [
        {
            "events": [
                {
                    "label": "Mix0",
                    "mixing_query": {"ch1": {"tip_angle": np.pi / 2, "phase": np.pi}},
                },
                {
                    "label": "Dur0",
                    "duration": 1,
                    "rotor_frequency": 0,
                    "rotor_angle": 0.5,
                    "transition_query": [{"ch1": {"P": [0], "D": [0]}}],
                },
                {
                    "label": "Spec0",
                    "fraction": 1,
                    "rotor_frequency": 20000,
                    "rotor_angle": 0.5,
                    "transition_query": [{"ch1": {"P": [1], "D": [0]}}],
                },
            ]
        },
        {
            "events": [
                {
                    "label": "Mix1",
                    "mixing_query": {"ch1": {"tip_angle": np.pi / 2, "phase": np.pi}},
                },
                {
                    "label": "Dur1",
                    "duration": 2,
                    "rotor_frequency": 0,
                    "rotor_angle": 0.5,
                    "transition_query": [{"ch1": {"P": [3], "D": [0]}}],
                },
                {
                    "label": "Spec1",
                    "fraction": 1,
                    "rotor_frequency": 0,
                    "rotor_angle": 0.5,
                    "transition_query": [{"ch1": {"P": [4], "D": [0]}}],
                },
            ]
        },
    ]

    method1 = Method(
        name="test-method-1",
        channels=["1H", "13C"],
        spectral_dimensions=all_defined_no_constant_spec_dims,
    )

    # Constant 'magnetic_flux_density' and 'rotor_frequency'
    method2 = Method(
        name="test-method-2",
        channels=["1H", "13C"],
        spectral_dimensions=const_mfd_rotor_ang_spec_dims,
    )

    assert method1.name == "test-method-1"
    assert len(method1.spectral_dimensions) == 2

    basic_summary_tests(method1)

    assert method2.name == "test-method-2"
    assert len(method2.spectral_dimensions) == 2

    args_summary_tests(method2)

    assert isinstance(method2.plot(), Figure)
