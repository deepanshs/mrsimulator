# -*- coding: utf-8 -*-
import pytest
from mrsimulator.methods import BlochDecaySpectrum


def test_01():
    error = "method requires exactly 1 channel"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        BlochDecaySpectrum(channels=["1H", "29Si"])
