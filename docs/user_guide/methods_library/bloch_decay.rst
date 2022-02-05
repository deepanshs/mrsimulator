Bloch Decay Spectrum
--------------------

The :py:class:`~mrsimulator.methods.BlochDecaySpectrum` class simulates the
Bloch decay spectrum.

.. testcode::

    from mrsimulator.methods import BlochDecaySpectrum

    method = BlochDecaySpectrum(
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

.. testcode::

  from mrsimulator.methods import BlochSpectrum

.. minigallery:: mrsimulator.methods.BlochDecaySpectrum mrsimulator.methods.BlochSpectrum
..     :add-heading: Examples using ``BlochDecaySpectrum``
..     :heading-level: "
