# -*- coding: utf-8 -*-
import json
from copy import deepcopy
from urllib.parse import urlparse

from astropy import units as u

from .__version__ import __version__
from ._utils_download_file import _download_file_from_url
from .simulator import _Isotopomers
from .simulator import _Spectrum
from .simulator import get_csdfpy_object
from .utils import __get_spin_attribute__


__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]
__version__ = __version__


def _fn_(x):
    return int("".join([i for i in x if i.isnumeric()]))


def _import_json(filename):
    res = urlparse(filename)
    if res[0] not in ["file", ""]:
        filename = _download_file_from_url(filename)
    with open(filename, "rb") as f:
        content = f.read()
        return json.loads(str(content, encoding="UTF-8"))


class Simulator:
    """
    The simulator class.
    """

    __slots__ = (
        "_isotopomers_c",
        "_isotopomers",
        "_spectrum",
        "_spectrum_c",
        "_isotope_list",
        "_allowed_isotopes",
        "_larmor_frequency",
        "_freq",
        "_amp",
    )

    def __init__(self, isotopomers=None, spectrum=None):
        self._isotopomers_c = []
        self._isotopomers = []
        self._spectrum = {}
        self._spectrum_c = {}
        self._isotope_list = []
        self._freq = [] * u.Unit("Hz")
        self._amp = []
        isotope_list = __get_spin_attribute__.keys()
        self._allowed_isotopes = list(
            set(
                [
                    isotope
                    for isotope in isotope_list
                    if __get_spin_attribute__[isotope]["spin"] == 0.5
                ]
            )
        )

        if isotopomers is not None:
            self.isotopomers = isotopomers

        if spectrum is not None:
            self.spectrum = spectrum

    @property
    def isotope_list(self):
        """
        Return a list of unique isotopes symbols from the list of
        isotopomers.
        """
        return self._isotope_list

    @property
    def isotopomers(self):
        """
        Return a list of :ref:`isotopomer` objects. The attribute can also be
        used to assign a list of valid isotopomers.
        """
        return self._isotopomers

    @isotopomers.setter
    def isotopomers(self, value):
        self._isotopomers_c = _Isotopomers(value)
        isotope_list = [
            site["isotope_symbol"]
            for isotopomer in self._isotopomers_c
            for site in isotopomer["sites"]
            if site["isotope_symbol"] in self._allowed_isotopes
        ]

        self._isotope_list = list(set(isotope_list))
        self._isotope_list.sort(key=_fn_)
        self._isotopomers = value

    @property
    def spectrum(self):
        """
        Return a :ref:`spectrum` object. The attribute can also be
        used to assign a valid spectrum object.
        """
        # return json.dumps(self._spectrum, ensure_ascii=True, indent=2)
        return self._spectrum

    @spectrum.setter
    def spectrum(self, value):
        self._spectrum_c = _Spectrum(**value["direct_dimension"])
        self._spectrum = value

    def run(self, method, data_object=False, **kwargs):
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


def run_test():
    from mrsimulator.methods import one_d_spectrum
    import matplotlib.pyplot as plt
    from . import examples

    s1 = Simulator()
    # test 1
    s1.isotopomers, s1.spectrum = examples.csa_static()
    freq, amp = s1.run(one_d_spectrum, verbose=1)

    ax = plt.subplots(1, 2, figsize=(6, 3))[1]
    ax[0].plot(freq, amp)
    # ax[0].plot(ob1.dimensions[0].coordinates,
    #            ob1.dependent_variables[0].components[0])
    # label_ = ob1.dimensions[0].axis_label
    label_ = f"frequency / {freq.unit}"
    ax[0].set_xlabel(label_)

    # test 2
    s1.isotopomers, s1.spectrum = examples.csa_mas()
    freq, amp = s1.run(one_d_spectrum, verbose=1)
    ax[1].plot(freq, amp)
    # ax[1].plot(ob2.dimensions[0].coordinates,
    #            ob2.dependent_variables[0].components[0])
    # label_ = ob1.dimensions[0].axis_label
    label_ = f"frequency / {freq.unit}"
    ax[1].set_xlabel(label_)
    plt.tight_layout()
    plt.show()
