# -*- coding: utf-8 -*-
from .unit import _ppm
from .unit import string_to_quantity
from .utils import __get_spin_attribute__

try:
    import csdfpy as cp
except ImportError:
    pass


__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


class _Dimensions:
    __slots__ = ()

    def __new__(
        number_of_points=1024,
        spectral_width="100 kHz",
        reference_offset="0 Hz",
    ):
        """Initialize."""
        dictionary = {}
        dictionary["number_of_points"] = int(number_of_points)

        spectral_width = string_to_quantity(spectral_width)
        if spectral_width.unit.physical_type != "frequency":
            raise Exception(
                ("A frequency value is required for the 'spectral_width'.")
            )
        dictionary["spectral_width"] = spectral_width.to("Hz").value

        reference_offset = string_to_quantity(reference_offset)
        if reference_offset.unit.physical_type != "frequency":
            raise Exception(
                ("A frequency value is required for the 'reference_offset'.")
            )
        dictionary["reference_offset"] = reference_offset.to("Hz").value
        return dictionary


class _Spectrum(_Dimensions):
    """Set up a virtual spin environment."""

    __slots__ = ()

    def __new__(
        self,
        magnetic_flux_density="9.4 T",
        rotor_frequency="0 kHz",
        rotor_angle="54.735 deg",
        rotor_phase="0 rad",
        nucleus="1H",
        *args,
        **kwargs,
    ):
        """Initialize"""
        dimension_dictionary = super(_Spectrum, self).__new__(*args, **kwargs)
        magnetic_flux_density = string_to_quantity(magnetic_flux_density)
        if magnetic_flux_density.unit.physical_type != "magnetic flux density":
            raise Exception(
                (
                    "A magnetic flux density quantity is required for "
                    "'magnetic_flux_density'."
                )
            )
        magnetic_flux_density = magnetic_flux_density.to("T").value

        rotor_frequency = string_to_quantity(rotor_frequency)
        if rotor_frequency.unit.physical_type != "frequency":
            raise Exception(
                ("A frequency quantity is required for 'rotor_frequency'.")
            )
        rotor_frequency = rotor_frequency.to("Hz").value

        rotor_angle = string_to_quantity(rotor_angle).to("rad").value
        rotor_phase = string_to_quantity(rotor_phase).to("rad").value

        dictionary = {
            "magnetic_flux_density": magnetic_flux_density,
            "rotor_frequency": rotor_frequency,
            "rotor_angle": rotor_angle,
            "rotor_phase": rotor_phase,
        }

        dictionary.update(dimension_dictionary)

        detect = get_proper_detector_nucleus(nucleus)
        try:
            spin_dictionary = __get_spin_attribute__[detect]
            spin_dictionary["isotope"] = detect
        except KeyError:
            raise Exception(f"Failed to simulates the {detect} spectrum.")

        dictionary.update(spin_dictionary)
        return dictionary