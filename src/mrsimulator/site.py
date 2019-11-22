# -*- coding: utf-8 -*-
"""Base Site class."""
from typing import ClassVar
from typing import Dict
from typing import Optional

from mrsimulator import AntisymmetricTensor
from mrsimulator import Parseable
from mrsimulator import SymmetricTensor
from mrsimulator.dimension import get_isotope_data
from pydantic import validator

__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"


class Site(Parseable):
    """
    Base Site class representing a nuclear isotope.

    Arguments:
        name: An optional string with a name or id of the site.
        isotope: An optional string expressed as atomic number followed by an
                isotope symbol, eg. ``13C``, ``17O``. The default value is ``1H``.
        isotropic_chemical_shift: An optional floating point number representing
                the isotropic chemical shift of the site in unit of ppm. The
                default value is ``0``.
        shielding_symmetric: An optional SymmetricTensor object or an equivalent
                python dict object representing the irreducible second-rank traceless
                symmetric part of the nuclear shielding tensor. The default value is
                ``None``.
        shielding_antisymmetric: An optional AntisymmetricTensor object or an
                equivalent python dict object representing the irreducible first-rank
                antisymmetric part of the nuclear shielding tensor. The default value
                is ``None``.
        quadrupolar: An optional SymmetricTensor object or an equivalent python dict
                object representing the irreducible second-rank traceless symmetric
                part of electric-field gradient tensor. The default value is ``None``.

    Example:
        Setting up Site objects.

            >>> site1 = Site(
            ...           isotope='13C',
            ...           isotropic_chemical_shift=20, # in ppm
            ...           shielding_symmetric={
            ...             "zeta": 10, # in ppm
            ...             "eta": 0.5
            ...           },
            ...           quadrupole={
            ...             "Cq": 5.1e6, # in Hz
            ...             "eta": 0.5
            ...           }
            ...         )

        Using SymmetricTensor objects.

        .. doctest::

            >>> site1 = Site(
            ...           isotope='13C',
            ...           isotropic_chemical_shift=20, # in ppm
            ...           shielding_symmetric=SymmetricTensor(zeta=10, eta=0.5),
            ...           quadrupole=SymmetricTensor(Cq=5.1e6, eta=0.5)
            ...         )
    """

    name: str = None
    isotope: str = "1H"
    isotropic_chemical_shift: Optional[float] = 0
    shielding_symmetric: Optional[SymmetricTensor] = None
    shielding_antisymmetric: Optional[AntisymmetricTensor] = None
    quadrupolar: Optional[SymmetricTensor] = None

    property_unit_types: ClassVar = {"isotropic_chemical_shift": "dimensionless"}
    property_default_units: ClassVar = {"isotropic_chemical_shift": "ppm"}
    property_units: Dict = {"isotropic_chemical_shift": "ppm"}

    @validator("quadrupolar")
    def spin_must_be_at_least_one(cls, v, values):
        isotope = values["isotope"]
        I = get_isotope_data(isotope)["spin"] / 2.0
        if I >= 1:
            if "zeta" in v.property_units:
                v.property_units.pop("zeta")
            return v
        message = f"with spin quantum number {I} does not allow quadrupolar tensor."
        raise ValueError(f"{isotope} {message}")

    @validator("shielding_symmetric")
    def shielding_symmetric_must_not_contain_Cq(cls, v, values):
        if "Cq" in v.property_units:
            v.property_units.pop("Cq")
        return v

    class Config:
        validate_assignment = True

    @property
    def spin(self):
        """Spin quantum number of the isotope."""
        return get_isotope_data(self.isotope)["spin"] / 2.0

    @property
    def natural_abundance(self):
        """Natural abundance of the isotope in units of %."""
        return get_isotope_data(self.isotope)["natural_abundance"]

    @property
    def gyromagnetic_ratio(self):
        """Reduced gyromagnetic ratio of the isotope in units of MHz/T."""
        return get_isotope_data(self.isotope)["gyromagnetic_ratio"]

    @property
    def quadrupole_moment(self):
        """Quadrupole moment of the isotope in units of eB (electron-Barn)."""
        return get_isotope_data(self.isotope)["quadrupole_moment"]

    @property
    def atomic_number(self):
        """Atomic number of the isotope."""
        return get_isotope_data(self.isotope)["atomic_number"]

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse the physical quantities from a dictionary representation of the Site
        object.

        Args:
            py_dict: A required Dict object.

        Returns:
            Site object

        Example:
            >>> site_dict = {
            ...    "isotope": "13C",
            ...    "isotropic_chemical_shift": "20 ppm",
            ...    "shielding_symmetric": {"zeta": "10 ppm", "eta":0.5}
            ... }
            >>> site1 = Site.parse_dict_with_units(site_dict)
        """
        prop_mapping = {
            "shielding_symmetric": SymmetricTensor,
            "shielding_antisymmetric": AntisymmetricTensor,
            "quadrupolar": SymmetricTensor,
        }

        for k, v in prop_mapping.items():
            if k in py_dict:
                py_dict[k] = v.parse_dict_with_units(py_dict[k])

        return super().parse_dict_with_units(py_dict)

    def to_freq_dict(self, B0):
        """
        Serialize the Site object to a JSON compliant python dictionary where the
        attribute values are numbers expressed in default units. The default unit
        for attributes with respective dimensionalities are:
        - frequency: `Hz`
        - angle: `rad`

        Args:
            B0: A required macroscopic magnetic flux density given in units of T.

        Return:
            Dict object

        Example:
            >>> pprint(site1.to_freq_dict(9.4))
            {'isotope': '13C',
             'isotropic_chemical_shift': -2013.1791999999998,
             'name': None,
             'quadrupolar': None,
             'shielding_antisymmetric': None,
             'shielding_symmetric': {'alpha': None,
                                     'beta': None,
                                     'eta': 0.5,
                                     'gamma': None,
                                     'zeta': -1006.5895999999999}}
        """
        temp_dict = self.dict()
        larmor_frequency = -self.gyromagnetic_ratio * B0  # in MHz
        for k in ["shielding_symmetric", "shielding_antisymmetric", "quadrupolar"]:
            if getattr(self, k):
                temp_dict[k] = getattr(self, k).to_freq_dict(larmor_frequency)
                if k == "shielding_symmetric":
                    temp_dict[k].pop("Cq")

        temp_dict["isotropic_chemical_shift"] *= larmor_frequency
        temp_dict.pop("property_units")
        return temp_dict

    def to_dict_with_units(self):
        """
        Serialize the Site object to a JSON compliant python dictionary with units.

        Returns:
            Dict object

        Example:
            >>> pprint(site1.to_dict_with_units())
            {'isotope': '13C',
             'isotropic_chemical_shift': '20.0 ppm',
             'shielding_symmetric': {'eta': 0.5, 'zeta': '10.0 ppm'}}

        """
        temp_dict = {k: v for k, v in self.dict().items() if v is not None}

        for k in ["shielding_symmetric", "shielding_antisymmetric", "quadrupolar"]:
            if getattr(self, k):
                temp_dict[k] = getattr(self, k).to_dict_with_units()

        temp_dict["isotropic_chemical_shift"] = (
            str(temp_dict["isotropic_chemical_shift"]) + " ppm"
        )
        temp_dict.pop("property_units")

        return temp_dict
