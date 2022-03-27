Multi-Quantum VAS
-----------------

The :py:class:`~mrsimulator.method.lib.ThreeQ_VAS`, :py:class:`~mrsimulator.method.lib.FiveQ_VAS`, and
:py:class:`~mrsimulator.method.lib.SevenQ_VAS` classes all simulate a multiple quantum VAS
spectrum. The spinning speed for all three methods is fixed at infinite speed. The spectrum
is also sheared such that the correlating dimensions are the isotropic dimension and the VAS dimension.

.. code-block:: python

    from mrsimulator.method.lib import ThreeQ_VAS
    from mrsimulator.method import SpectralDimension

    method = ThreeQ_VAS(
        channels=["87Rb"],
        magnetic_flux_density=7,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=128,
                spectral_width=3e3,  # in Hz
                reference_offset=-2e3,  # in Hz
                label="Isotropic dimension",
            ),
            SpectralDimension(
                count=512,
                spectral_width=1e4,  # in Hz
                reference_offset=-5e3,  # in Hz
                label="MAS dimension",
            ),
        ],
    )

The other methods, representing 5 and 7 quantum transitions, can be imported as follows:

.. code-block:: python

    from mrsimulator.method.lib import FiveQ_VAS
    from mrsimulator.method.lib import SevenQ_VAS

.. minigallery:: mrsimulator.method.lib.ThreeQ_VAS mrsimulator.method.lib.FiveQ_VAS mrsimulator.method.lib.SevenQ_VAS
    :add-heading: Examples using the Multi-quantum VAS method.
    :heading-level: "
