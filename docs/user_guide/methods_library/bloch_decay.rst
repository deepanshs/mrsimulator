Bloch Decay Spectrum
--------------------

The :py:class:`~mrsimulator.methods.BlochDecaySpectrum` class simulates the
Bloch decay spectrum.

.. code-block:: python

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

.. minigallery:: mrsimulator.methods.BlochDecaySpectrum
    :add-heading: Examples using ``BlochDecaySpectrum``
    :heading-level: "
