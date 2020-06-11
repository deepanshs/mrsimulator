# -*- coding: utf-8 -*-
"""Base Site class."""
from typing import ClassVar
from typing import Dict
from typing import Optional

from pydantic import validator

from .isotope import Isotope
from .parseable import Parseable
from .tensors import AntisymmetricTensor
from .tensors import SymmetricTensor

__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"


class Site(Parseable):
    """
    Base class representing a single-site nuclear spin interaction tensor parameters.

    Attributes:
        name: An optional string with a name or id of the site. The default value is
            None.
        label: An optional label for the site. The default value is None.
        description: An optional description of the site. The default value is None.
        isotope: An optional string expressed as atomic number followed by an
            isotope symbol, eg. `13C`, `17O`. The default value is `1H`.
        isotropic_chemical_shift: An optional floating point number representing
            the isotropic chemical shift of the site in unit of ppm. The default value
            is 0.
        shielding_symmetric: An optional SymmetricTensor object or an equivalent python
            dict object representing the irreducible second-rank traceless symmetric
            part of the nuclear shielding tensor. The default value is None.
        shielding_antisymmetric: An optional AntisymmetricTensor object or an
            equivalent python dict object representing the irreducible first-rank
            antisymmetric part of the nuclear shielding tensor. The default value is
            None.
        quadrupolar: An optional SymmetricTensor object or an equivalent python dict
            object representing the irreducible second-rank traceless symmetric part of
            the electric-field gradient tensor. The default value is None.

    Example
    -------

    Setting up Site objects.

    >>> site1 = Site(
    ...     isotope='13C',
    ...     isotropic_chemical_shift=20, # in ppm
    ...     shielding_symmetric={
    ...         "zeta": 10, # in ppm
    ...         "eta": 0.5
    ...     },
    ...     quadrupole={
    ...         "Cq": 5.1e6, # in Hz
    ...         "eta": 0.5
    ...     }
    ... )

    Using SymmetricTensor objects.

    >>> site1 = Site(
    ...     isotope='13C',
    ...     isotropic_chemical_shift=20, # in ppm
    ...     shielding_symmetric=SymmetricTensor(zeta=10, eta=0.5),
    ...     quadrupole=SymmetricTensor(Cq=5.1e6, eta=0.5)
    ... )
    """

    name: str = None
    label: str = None
    description: str = None
    isotope: Optional[str] = "1H"
    isotropic_chemical_shift: Optional[float] = 0
    shielding_symmetric: Optional[SymmetricTensor] = None
    shielding_antisymmetric: Optional[AntisymmetricTensor] = None
    quadrupolar: Optional[SymmetricTensor] = None

    property_unit_types: ClassVar = {"isotropic_chemical_shift": "dimensionless"}
    property_default_units: ClassVar = {"isotropic_chemical_shift": "ppm"}
    property_units: Dict = {"isotropic_chemical_shift": "ppm"}

    @validator("quadrupolar")
    def spin_must_be_at_least_one(cls, v, values):
        if v is None:
            return v
        isotope = values["isotope"]
        isotope = Isotope(symbol=isotope) if isinstance(isotope, str) else isotope
        spin_I = isotope.spin
        if spin_I >= 1:
            if "zeta" in v.property_units:
                v.property_units.pop("zeta")
            return v
        message = (
            f"with spin quantum number {spin_I} does not allow quadrupolar tensor."
        )
        raise ValueError(f"{isotope} {message}")

    @validator("shielding_symmetric")
    def shielding_symmetric_must_not_contain_Cq(cls, v, values):
        if "Cq" in v.property_units:
            v.property_units.pop("Cq")
        return v

    @validator("isotope", always=True)
    def validate_isotope(cls, v, *, values, **kwargs):
        return Isotope(symbol=v)

    class Config:
        validate_assignment = True

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse the physical quantities from a dictionary representation of the Site
        object.

        Args:
            dict py_dict: A required Dict object.

        Returns:
            Site object

        Example
        -------

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
            float B0: A required macroscopic magnetic flux density given in units of T.

        Return:
            Dict object

        Example
        -------

        >>> pprint(site1.to_freq_dict(9.4))
        {'description': None,
         'isotope': '13C',
         'isotropic_chemical_shift': -2013.1791999999998,
         'label': None,
         'name': None,
         'quadrupolar': None,
         'shielding_antisymmetric': None,
         'shielding_symmetric': {'alpha': None,
                                 'beta': None,
                                 'eta': 0.5,
                                 'gamma': None,
                                 'zeta': -1006.5895999999999}}
        """
        temp_dict = self.dict(exclude={"isotope"})
        temp_dict["isotope"] = self.isotope.symbol
        larmor_frequency = -self.isotope.gyromagnetic_ratio * B0  # in MHz
        for k in ["shielding_symmetric", "shielding_antisymmetric", "quadrupolar"]:
            if getattr(self, k):
                temp_dict[k] = getattr(self, k).to_freq_dict(larmor_frequency)
                if k == "shielding_symmetric":
                    temp_dict[k].pop("Cq")

        temp_dict["isotropic_chemical_shift"] *= larmor_frequency
        temp_dict.pop("property_units")
        return temp_dict
