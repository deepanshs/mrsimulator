SSB2D
-----

The :py:class:`~mrsimulator.method.lib.SSB2D` class simulates a sheared 2D finite
to infinite speed MAS correlation spectrum. The spinning speed for the second spectral
dimension is fixed at infinite spinning speed

.. code-block:: python

    from mrsimulator.method.lib import SSB2D
    from mrsimulator.method import SpectralDimension

    method = SSB2D(
        channels=["13C"],
        magnetic_flux_density=7,  # in T
        rotor_frequency=1500,  # in Hz
        spectral_dimensions=[
            SpectralDimension(
                count=16,
                spectral_width=16 * 1500,  # in Hz (= count * rotor_frequency)
                reference_offset=-5e3,  # in Hz
                label="Sideband dimension",
            ),
            SpectralDimension(
                count=512,
                spectral_width=1e4,  # in Hz
                reference_offset=-4e3,  # in Hz
                label="Isotropic dimension",
            ),
        ],
    )

**Examples**

.. minigallery:: mrsimulator.method.lib.SSB2D
