# -*- coding: utf-8 -*-
from copy import deepcopy

from astropy import units as u

from mrsimulator.utils import _download_file_from_url, _fn_, _import_json
from mrsimulator import Isotopomer, Spectrum
from mrsimulator.spectrum import ISOTOPE_DATA

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


class Simulator:
    """
    The simulator class.
    """

    def __init__(self, isotopomers=None):
        self.isotopomers = isotopomers or []
        self.isotope_list = []

    @staticmethod
    def allowed_isotopes():
        """
        Returns a list of all valid isotopes for this simulator
        """
        return list(
            {
                isotope
                for isotope in isotope_list
                if ISOTOPE_DATA[isotope]["spin"] == 1
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
                site.nucleus
                for isotopomer in isotopomers
                for site in isotopomer.sites
            }
        )

    @property
    def valid_isotope_list(self):
        """
        Returns a list of unique and valid isotope symbols from the list of isotopomers
        """
        return list(
            {
                site.nucleus
                for isotopomer in isotopomers
                for site in isotopomer.sites
                if site.nucleus in self.allowed_isotopes
            }
        )

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

        isotopomers = self._isotopomers_c
        if isotopomers is []:
            raise Exception("Isotopomers are required for simulation.")
        spectrum = self._spectrum_c
        if spectrum is {}:
            raise Exception(
                ("Cannot simulate without the spectrum information.")
            )
        (
            self._freq,
            self._amp,
            self._larmor_frequency,
            list_index_isotopomer,
        ) = method(spectrum=spectrum, isotopomers=isotopomers, **kwargs)
        """The frequency is in the units of Hz."""
        self._freq *= u.Unit("Hz")
        """The larmor_frequency is in the units of MHz."""
        self._larmor_frequency *= u.Unit("MHz")

        isotopo_ = [isotopomers[i] for i in list_index_isotopomer]
        application = {"isotopomers": str(isotopo_), "spectrum": spectrum}
        data_object = False
        if data_object:
            return get_csdfpy_object(
                self._freq, self._larmor_frequency, self._amp, application
            )
        else:
            return deepcopy(self._freq), deepcopy(self._amp)

    def load_isotopomers(self, filename):
        """
        Load a JSON serialized isotopomers file.

        See an
        `example <https://raw.githubusercontent.com/DeepanshS/mrsimulator-test
        /master/isotopomers_ppm.json>`_
        of JSON serialized isotopomers file. For details, refer to the
        :ref:`load_isotopomers` section.
        """
        contents = _import_json(filename)
        self.isotopomers = contents["isotopomers"]
