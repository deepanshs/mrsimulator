Satellite-Transition MAS
------------------------

The :py:class:`~mrsimulator.methods.ST1_VAS` and :py:class:`~mrsimulator.methods.ST2_VAS` classes
simulate a sheared and scaled satellite and central transition correlation spectrum. The spinning
speed for these methods is fixed at infinite speed.

.. testcode::

    from mrsimulator.methods import ST1_VAS

    method = ST1_VAS(
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            dict(
                count=128,
                spectral_width=1e3,  # in Hz
                reference_offset=-5e3,  # in Hz
                label="Isotropic dimension",
            ),
            dict(
                count=256,
                spectral_width=1e4,  # in Hz
                reference_offset=-3e3,  # in Hz
                label="MAS dimension",
            ),
        ],
    )

.. minigallery:: mrsimulator.methods.ST1_VAS mrsimulator.methods.ST2_VAS
    :add-heading: Examples using ``ST1_VAS`` and ``ST2_VAS``
    :heading-level: "
