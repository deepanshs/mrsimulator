import numpy as np
import pytest
from mrsimulator.method import SpectralDimension
from mrsimulator.method.lib.base import BaseNamedMethod1D
from mrsimulator.method.lib.base import BaseNamedMethod2D
from mrsimulator.utils.error import ImmutableEventError
from mrsimulator.utils.error import NamedMethodError


def test_BaseNamedMethod1D_spectral_dimension_count():
    e = "Method requires exactly 1 spectral dimensions, given 2."
    # test for SpectralDimension passed as a dictionary
    with pytest.raises(ValueError, match=f".*{e}.*"):
        BaseNamedMethod1D(spectral_dimensions=[{}, {}])
    # test for SpectralDimension passed as object
    with pytest.raises(ValueError, match=f".*{e}.*"):
        BaseNamedMethod1D(
            spectral_dimensions=[
                SpectralDimension(),
                SpectralDimension(),
            ]
        )


def test_BaseNamedMethod1D_setting_events():
    e = "Event objects are immutable for BaseNamedMethod1D class."
    with pytest.raises(ImmutableEventError, match=f".*{e}.*"):
        BaseNamedMethod1D(spectral_dimensions=[{"events": [{"rotor_angle": 1.32}]}])


def test_BaseNamedMethod1D_setting_name():
    e = "`name=test != classname=BaseNamedMethod1D`."
    with pytest.raises(NamedMethodError, match=f".*{e}.*"):
        BaseNamedMethod1D(channels=["1H"], name="test", spectral_dimensions=[{}])


def test_BaseNamedMethod2D_spectral_dimension_count():
    e = "Method requires exactly 2 spectral dimensions, given 1."
    # test for SpectralDimension passed as a dictionary
    with pytest.raises(ValueError, match=f".*{e}.*"):
        BaseNamedMethod2D(channels=["1H"], spectral_dimensions=[{}])
    # test for SpectralDimension passed as object
    with pytest.raises(ValueError, match=f".*{e}.*"):
        BaseNamedMethod2D(channels=["1H"], spectral_dimensions=[SpectralDimension()])


def test_BaseNamedMethod2D_setting_events():
    e = "Event objects are immutable for BaseNamedMethod2D class."
    with pytest.raises(ImmutableEventError, match=f".*{e}.*"):
        BaseNamedMethod2D(
            channels=["1H"],
            spectral_dimensions=[{"events": [{"rotor_angle": 1.32}]}, {}],
        )


def test_BaseNamedMethod2D_setting_name():
    e = "`name=test != classname=BaseNamedMethod2D`."
    with pytest.raises(NamedMethodError, match=f".*{e}.*"):
        BaseNamedMethod2D(channels=["1H"], name="test", spectral_dimensions=[{}, {}])


def test_BaseNamedMethod_wrong_number_of_channels():
    e = r"BaseNamedMethod1D only supports 1 channel\(s\). Got 2 channels"
    with pytest.raises(ValueError, match=e):
        BaseNamedMethod1D(channels=["1H", "13C"], name="BaseNamedMethod1D")

    e = r"BaseNamedMethod2D only supports 1 channel\(s\). Got 2 channels"
    with pytest.raises(ValueError, match=e):
        BaseNamedMethod2D(channels=["1H", "13C"], name="BaseNamedMethod2D")


def test_BaseNamedMethod2D_set_inf_rotor_freq():
    mtd = BaseNamedMethod2D(
        channels=["17O"], name="BaseNamedMethod2D", rotor_frequency=np.inf
    )
    assert mtd.rotor_frequency == 1e12
