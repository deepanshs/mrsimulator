Satellite-Transition VAS
------------------------

The :py:class:`~mrsimulator.method.lib.ST1_VAS` and :py:class:`~mrsimulator.method.lib.ST2_VAS` classes
simulate a sheared and scaled satellite and central transition correlation spectrum. The spinning
speed for these methods is fixed at infinite speed.

.. code-block:: python

    from mrsimulator.method.lib import ST1_VAS
    from mrsimulator.method import SpectralDimension

    method = ST1_VAS(
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=128,
                spectral_width=1e3,  # in Hz
                reference_offset=-5e3,  # in Hz
                label="Isotropic dimension",
            ),
            SpectralDimension(
                count=256,
                spectral_width=1e4,  # in Hz
                reference_offset=-3e3,  # in Hz
                label="MAS dimension",
            ),
        ],
    )

.. minigallery:: mrsimulator.method.lib.ST1_VAS mrsimulator.method.lib.ST2_VAS
