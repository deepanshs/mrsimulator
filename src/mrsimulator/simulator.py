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
    def allowed_isotopes():
        """
        Returns a list of all valid isotopes for this simulator
        """
        return list({isotope for isotope, data in ISOTOPE_DATA.items()})

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

    @property
    def unique_isotopes(self):
        """
        Returns a list of unique and valid isotope symbols from the list of isotopomers
        """
        return list(
            {
                site.isotope
                for isotopomer in self.isotopomers
                for site in isotopomer.sites
                if site.isotope in self.allowed_isotopes()
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
        spectrum = self.spectrum.dict()

        (freq, amp, larmor_frequency, list_index_isotopomer) = one_d_spectrum(
            spectrum=spectrum, isotopomers=isotopomers, **kwargs
        )
        """The frequency is in the units of Hz."""
        freq *= u.Unit("Hz")
        """The larmor_frequency is in the units of MHz."""
        # larmor_frequency *= u.Unit("MHz")

        # isotopo_ = [self.isotopomers[i] for i in list_index_isotopomer]
        return freq, amp
