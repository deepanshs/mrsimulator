# -*- coding: utf-8 -*-
import warnings

warnings.warn(
    (
        "Importing library methods from `mrsimulator.methods` is deprectated and will "
        "be removed in the next version. Please import library methods from the "
        "`mrsimulator.method.lib` module."
    ),
    Warning,
)

from mrsimulator.method.lib import BlochDecayCentralTransitionSpectrum  # noqa:F401
from mrsimulator.method.lib import BlochDecayCTSpectrum  # noqa:F401
from mrsimulator.method.lib import BlochDecaySpectrum  # noqa:F401
from mrsimulator.method.lib import FiveQ_VAS  # noqa:F401
from mrsimulator.method.lib import SevenQ_VAS  # noqa:F401
from mrsimulator.method.lib import ThreeQ_VAS  # noqa:F401
from mrsimulator.method.lib import SSB2D  # noqa:F401
from mrsimulator.method.lib import ST1_VAS  # noqa:F401
from mrsimulator.method.lib import ST2_VAS  # noqa:F401

# This file is intended only for raising DepreciationWarning on importing library
# methods from mrsimulator.methods path
# REMOVE THIS FILE UPON v0.8 RELEASE
