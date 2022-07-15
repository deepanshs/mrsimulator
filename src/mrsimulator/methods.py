import warnings

warnings.warn(
    (
        "Importing library methods from `mrsimulator.methods` is deprecated and will "
        "be removed in the next version. Please import library methods from the "
        "`mrsimulator.method.lib` module."
    ),
    Warning,
)

from mrsimulator.method.lib import *  # noqa:F401,F403

# This file is intended only for raising DepreciationWarning on importing library
# methods from mrsimulator.methods path
# REMOVE THIS FILE UPON v0.8 RELEASE
