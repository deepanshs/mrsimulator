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
        "isotropic_chemical_shift": ["dimensionless", "frequency"]
    }

    property_default_units: ClassVar = {
        "isotropic_chemical_shift": ["ppm", "Hz"]
    }

    @classmethod
    def parse_json_with_units(cls, json_dict):

        prop_mapping = {
            "shielding_symmetric": SymmetricTensor,
            "shielding_antisymmetric": AntisymmetricTensor,
            "quadrupolar": SymmetricTensor,
        }

        for k, v in prop_mapping.items():
            if k in json_dict:
                json_dict[k] = v.parse_json_with_units(json_dict[k])

        return super().parse_json_with_units(json_dict)

    def to_freq_dict(self, larmor_frequency):
        """
        Enforces units of Hz by multiplying any ppm values by the Larmor frequency in MHz
        MHz*ppm -> Hz
        """
        temp_dict = self.dict()

        for k in [
            "shielding_symmetric",
            "shielding_antisymmetric",
            "quadrupolar",
        ]:
            if getattr(self,k):
                temp_dict[k] = getattr(self,k).to_freq_dict(larmor_frequency)

        if self.property_units["isotropic_chemical_shift"] is "ppm":
            temp_dict["isotropic_chemical_shift"] *= larmor_frequency

        return temp_dict
