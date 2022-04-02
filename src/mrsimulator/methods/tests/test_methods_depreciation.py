# -*- coding: utf-8 -*-
import pytest


def test_import():
    warning = (
        "Importing library methods from `mrsimulator.methods` is deprectated and will "
        "be removed in the next version. Please import library methods from the "
        "`mrsimulator.method.lib` module."
    )
    with pytest.warns(Warning, match=f".*{warning}.*"):
        import mrsimulator.methods  # noqa:F401

        # from mrsimulator.methods import BlochDecaySpectrum
