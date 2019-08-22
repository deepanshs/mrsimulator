# -*- coding: utf-8 -*-
# from pydash import has, get
from astropy import units as u
from mrsimulator import Isotopomer, Spectrum
from mrsimulator.spectrum import ISOTOPE_DATA
from mrsimulator.methods import one_d_spectrum
from mrsimulator.importer import import_json


__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


class Simulator:
    """
    The simulator class.
    """

    def __init__(self, isotopomers=[], spectrum={}):
        self.isotopomers = isotopomers
        self.spectrum = spectrum

    @staticmethod
    def allowed_isotopes(spin=None):
        """
        Returns a list of all valid isotopes for this simulator
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
    def all_isotopes(self):
        """
        Return a list of unique isotopes symbols from the list of
        isotopomers.
        """
        return list(
            {
                site.isotope
                for isotopomer in self.isotopomers
                for site in isotopomer.sites
            }
        )

    def isotope_list(self, spin=None):
        """
        Returns a list of unique and valid isotope symbols from the list of isotopomers
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
        Load a JSON serialized isotopomers file.

        See an
        `example <https://raw.githubusercontent.com/DeepanshS/mrsimulator-test
        /master/isotopomers_ppm.json>`_
        of JSON serialized isotopomers file. For details, refer to the
        :ref:`load_isotopomers` section.
        """
        contents = import_json(filename)
        json_data = contents["isotopomers"]
        self.isotopomers = [Isotopomer.parse_json_with_units(obj) for obj in json_data]

    def run(self, method, **kwargs):
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
            isotopomer.to_freq_dict(self.spectrum.larmor_frequency)
            for isotopomer in self.isotopomers
        ]
        spectrum = self.spectrum

        freq, amp = one_d_spectrum(spectrum=spectrum, isotopomers=isotopomers, **kwargs)
        """The frequency is in the units of Hz."""
        freq *= u.Unit("Hz")
        return freq, amp
