Bloch Decay Central Transition
------------------------------

The :py:class:`~mrsimulator.method.lib.BlochDecayCTSpectrum` class simulates the
Bloch decay central transition selective spectrum.

.. code-block:: python

    from mrsimulator.method.lib import BlochDecayCTSpectrum
    from mrsimulator.method import SpectralDimension

    method = BlochDecayCTSpectrum(
        channels=["87Rb"],
        rotor_frequency=12500,  # in Hz
        rotor_angle=54.735 * 3.14159 / 180,  # in rad
        magnetic_flux_density=9.4,  # in tesla
        spectral_dimensions=[
            SpectralDimension(
                count=1024,
                spectral_width=25e3,  # in Hz
                reference_offset=-4e3,  # in Hz
            )
        ],
    )

**Examples**

.. minigallery:: mrsimulator.method.lib.BlochDecayCTSpectrum
