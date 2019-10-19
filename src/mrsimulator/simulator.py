# -*- coding: utf-8 -*-
# from pydash import has, get
from astropy import units as u
from mrsimulator import Isotopomer
from mrsimulator.spectrum import ISOTOPE_DATA
from mrsimulator.methods import one_d_spectrum
from mrsimulator.importer import import_json
import csdmpy as cp

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


class Simulator:
    """
    The simulator class.

    .. rubric:: Attributes Documentation

    Attributes:
        isotopomers: List of Isotopomer objects.
        spectrum: List of dimension objects.
        isotope: List of all unique isotopes defined in the list of isotopomers.
            This also includes NMR inactive isotopes.
    """

    def __init__(self, isotopomers=[], spectrum=[]):
        self.isotopomers = isotopomers
        self.spectrum = spectrum
        self.simulated_data = []

    @staticmethod
    def allowed_isotopes(spin=None):
        """
        List of NMR active isotopes allowed in ``mrsimulator``.

        Args:
            spin: (optional) The spin quantum number. Valid input are multiples of 0.5.

        Returns:
            A list of all isotopes with the give spin quantum number allowed in
            mrsimulator. If the spin is unspecified or None, a list of all
            allowed isotopes is returned instead.
        """
        if spin is None:
            return list({isotope for isotope, data in ISOTOPE_DATA.items()})
        return list(
            {
                isotope
                for isotope, data in ISOTOPE_DATA.items()
                if data["spin"] == int(2 * spin)
            }
        )

    @property
    def isotopes(self):
        return list(
            {
                site.isotope
                for isotopomer in self.isotopomers
                for site in isotopomer.sites
            }
        )

    def get_isotopes(self, spin=None):
        """
        List of unique isotopes defined in list of isotopomers.

        Args:
            spin: (optional) The spin quantum number. Valid input are multiples of 0.5.

        Returns:
            A list of unique isotopes from the list of isotopomers corresponding
            to given value of `spin`. If the spin is unspecified or None, a list of all
            unique isotopes is returned instead.
        """
        return list(
            {
                site.isotope
                for isotopomer in self.isotopomers
                for site in isotopomer.sites
                if site.isotope in self.allowed_isotopes(spin)
            }
        )

    def load_isotopomers(self, filename):
        """
        Load a list of isotopomers from JSON serialized isotopomers file.

        See an
        `example <https://raw.githubusercontent.com/DeepanshS/mrsimulator-test
        /master/isotopomers_ppm.json>`_
        of JSON serialized isotopomers file. For details, refer to the
        :ref:`load_isotopomers` section.

        Args:
            `filename`: A local or remote address to the JSON serialized isotopomers
                        file.
        """
        contents = import_json(filename)
        json_data = contents["isotopomers"]
        self.isotopomers = [Isotopomer.parse_dict_with_units(obj) for obj in json_data]

    def run(self, method, **kwargs):
        """Simulate the lineshape.

        Args:
            method: The methods used in simulation of the line-shapes.
        """
        return self.one_d_spectrum(**kwargs)

    def one_d_spectrum(self, **kwargs):
        """
        Simulate the spectrum using the specified method. The keyword argument
        are the arguments of the specified `method`.

        :ivar method: The method used in computing the linshape.
        :ivar data_object: If true, returns a `csdm` data object. If false,
            returns a tuple of frequency array and the
            corresponding amplitude array. The amplitude is a
            `numpy <https://docs.scipy.org/doc/numpy/reference/generated
            /numpy.array.html>`_
            array, whereas, the frequency is a
            `Quantity <http://docs.astropy.org/en/stable/units/quantity.html>`_
            array. The default value is False.
        """

        isotopomers = [
            isotopomer.to_freq_dict(self.spectrum[0].larmor_frequency)
            for isotopomer in self.isotopomers
        ]
        spectrum = self.spectrum

        freq, amp = one_d_spectrum(
            spectrum=spectrum[0], isotopomers=isotopomers, **kwargs
        )
        self.simulated_data = amp

        """The frequency is in the units of Hz."""
        freq *= u.Unit("ppm")
        return freq, amp

    def as_csdm_object(self):
        new = cp.new()
        for spectrum in self.spectrum:
            larmor_frequency = spectrum.larmor_frequency
            dimension = {
                "type": "linear",
                "count": spectrum.number_of_points,
                "increment": "{0} Hz".format(
                    spectrum.spectral_width / spectrum.number_of_points
                ),
                "coordinates_offset": f"{spectrum.reference_offset} Hz",
                "origin_offset": f"{larmor_frequency} MHz",
                "complex_fft": True,
            }
            new.add_dimension(dimension)

        if isinstance(self.simulated_data, list):
            for index, datum in enumerate(self.simulated_data):
                if datum != []:
                    dependent_variable = {
                        "type": "internal",
                        "quantity_type": "scalar",
                        "numeric_type": "float64",
                        "components": [datum],
                    }

                    name = self.isotopomers[index].name
                    if name != "":
                        dependent_variable.update({"name": name})

                    description = self.isotopomers[index].description
                    if description != "":
                        dependent_variable.update({"description": description})

                    dependent_variable["application"] = {
                        "com.github.DeepanshS.mrsimulator": {
                            "isotopomers": [
                                self.isotopomers[index].to_dict_with_units()
                            ]
                        }
                    }

                    new.add_dependent_variable(dependent_variable)
                    new.dependent_variables[-1].encoding = "base64"

        else:
            dependent_variable = {
                "type": "internal",
                "quantity_type": "scalar",
                "numeric_type": "float64",
                "components": [self.simulated_data],
            }
            new.add_dependent_variable(dependent_variable)
            new.dependent_variables[-1].encoding = "base64"

        return new
