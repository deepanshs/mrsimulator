Bloch Decay Central Transition
------------------------------

The :py:class:`~mrsimulator.methods.BlochDecayCTSpectrum` class simulates the
Bloch decay central transition selective spectrum.

.. code-block:: python

    from mrsimulator.methods import BlochDecayCTSpectrum

    method = BlochDecayCTSpectrum(
        channels=["1H"],
        rotor_frequency=12500,  # in Hz
        rotor_angle=54.735 * 3.14159 / 180,  # in rad
        magnetic_flux_density=9.4,  # in tesla
        spectral_dimensions=[
            dict(
                count=1024,
                spectral_width=25e3,  # in Hz
                reference_offset=-4e3,  # in Hz
            )
        ],
    )

The method may be imported using the following alias classes:

.. code-block:: python

    from mrsimulator.methods import BlochDecayCentralTransitionSpectrum
    from mrsimulator.methods import BlochCTSpectrum

.. minigallery:: mrsimulator.methods.BlochDecayCTSpectrum mrsimulator.methods.BlochDecayCentralTransitionSpectrum mrsimulator.methods.BlochCTSpectrum
    :add-heading:
    :heading-level: "
