# -*- coding: utf-8 -*-
from typing import ClassVar, Optional
from mrsimulator import Parseable, SymmetricTensor, AntisymmetricTensor

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


class Site(Parseable):
    """
    Base Site class representing a nuclear spin.

    .. rubric:: Attributes Documentation

    Attributes:
        isotope: A required string expressed as atomic number followed by an
                isotope symbol, eg. "13C", "17O". The default is "1H".
        isotropic_chemical_shift: An optional floating point number representing
                the isotropic chemical shift of the site in units of ppm. The
                default value is 0.
        shielding_symmetric: An optional SymmetricTensor object representing the
                traceless symmetric second-rank nuclear shielding tensor. The
                default value is None.
        shielding_antisymmetric: An optional AntisymmetricTensor object
                representing the antisymmetric first-rank nuclear shielding tensor.
                The default value is None.
        quadrupolar: An optional SymmetricTensor object representing the traceless
                symmetric second-rank electric-field gradient tensor. The
                default value is None.
    """

    isotope: str = "1H"
    isotropic_chemical_shift: Optional[float] = 0
    shielding_symmetric: Optional[SymmetricTensor] = None
    shielding_antisymmetric: Optional[AntisymmetricTensor] = None
    quadrupolar: Optional[SymmetricTensor] = None

    property_unit_types: ClassVar = {"isotropic_chemical_shift": "dimensionless"}

    property_default_units: ClassVar = {"isotropic_chemical_shift": "ppm"}

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse the physical quantities of a Site object when expressed as a
        python dictionary.

        Args:
            py_dict: Python dictionary representation of an isotopomers with
                        physical quantities.
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

    def to_freq_dict(self, larmor_frequency):
        """
        Enforces units of Hz by multiplying any ppm values by the Larmor frequency in
        MHz, MHz*ppm -> Hz
        """
        temp_dict = self.dict()

        for k in ["shielding_symmetric", "shielding_antisymmetric", "quadrupolar"]:
            if getattr(self, k):
                temp_dict[k] = getattr(self, k).to_freq_dict(larmor_frequency)

        temp_dict["isotropic_chemical_shift"] *= larmor_frequency

        return temp_dict

    def to_dict_with_units(self):
        """
        Serialize the Site object to a JSON compliant python dictionary with units.
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
