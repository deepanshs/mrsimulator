# -*- coding: utf-8 -*-
from copy import deepcopy

from astropy import units as u

from mrsimulator.utils import _download_file_from_url, _fn_, _import_json
from mrsimulator import Isotopomer, Spectrum
from mrsimulator.spectrum import ISOTOPE_DATA
from mrsimulator.methods import one_d_spectrum

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


class Simulator:
    """
    The simulator class.
    """

    def __init__(self, isotopomers, spectrum):
        self.isotopomers = isotopomers
        self.spectrum = spectrum

    @staticmethod
    def allowed_isotopes():
        """
        Returns a list of all valid isotopes for this simulator
        """
        return list({isotope for isotope, data in ISOTOPE_DATA.items() if data["spin"] == 1})

    @property
    def all_isotopes(self):
        """
        Return a list of unique isotopes symbols from the list of
        isotopomers.
        """
        return list({site.isotope for isotopomer in self.isotopomers for site in isotopomer.sites})

    @property
    def valid_isotope_list(self):
        """
        Returns a list of unique and valid isotope symbols from the list of isotopomers
        """
        return list({
            site.isotope
            for isotopomer in self.isotopomers for site in isotopomer.sites if site.isotope in self.allowed_isotopes()
        })

    @property
    def one_d_spectrum(self):
        """
        Get's a 1D spectrum for this
        """
        return self.run(one_d_spectrum)
    

    def run(self, method, **kwargs):
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

        :returns: A `csdm` object if `data_object` is True, else a tuple of
            frequency and amplitude. For details, refer to the
            description of `data_object`.
        """

        isotopomers = [isotopomer.dict() for isotopomer in self.isotopomers]
        spectrum = self.spectrum.dict()
        (
            freq,
            amp,
            larmor_frequency,
            list_index_isotopomer,
        ) = method(
            spectrum=spectrum, isotopomers=isotopomers, **kwargs)
        """The frequency is in the units of Hz."""
        freq *= u.Unit("Hz")
        """The larmor_frequency is in the units of MHz."""
        larmor_frequency *= u.Unit("MHz")

        isotopo_ = [self.isotopomers[i] for i in list_index_isotopomer]
        return freq, amp
