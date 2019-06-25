from typing import ClassVar, Optional
from mrsimulator import Parseable, SymmetricTensor, AntisymmetricTensor

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


class Site(Parseable):

    isotope: str = "1H"
    isotropic_chemical_shift: Optional[float]
    shielding_symmetric: Optional[SymmetricTensor]
    shielding_antisymmetric: Optional[AntisymmetricTensor]
    quadrupolar: Optional[SymmetricTensor]

    property_unit_types: ClassVar = {
        "isotropic_chemical_shift": ["dimensionless","frequency"],
        "anisotropy": ["dimensionless","frequency"],
        "alpha": "angle",
        "beta": "angle",
        "gamma": "angle",
    }

    property_default_units: ClassVar = {
        "isotropic_chemical_shift": ["ppm","Hz"],
    }

    @classmethod
    def parse_json_with_units(cls, json_dict):



        return super().parse_json_with_units(json_dict)
