# -*- coding: utf-8 -*-
from enum import Enum

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

freq_list_all = [
    "Shielding1_0",
    "Shielding1_2",
    "Quad1_2",
    "Quad2_0",
    "Quad2_2",
    "Quad2_4",
]


class FrequencyEnum(str, Enum):
    """Enumeration for selecting specific frequency contributions. The enumerations
    are:

    Attributes
    ----------

    Shielding1_0:
        Selects first-order and zeroth-rank nuclear shielding frequency contributions.

    Shielding1_2:
        Selects first-order and second-rank nuclear shielding frequency contributions.

    Quad1_2:
        Selects first-order and second-rank quadrupolar frequency contributions.

    Quad2_0:
        Selects second-order and zeroth-rank quadrupolar frequency contributions.

    Quad2_2:
        Selects second-order and second-rank quadrupolar frequency contributions.

    Quad2_4:
        Selects second-order and fourth-rank quadrupolar frequency contributions.
    """

    # Shielding1:
    #     Selects first-order nuclear shielding frequency contributions. Equivalent to
    #     ["Shielding1_0", "Shielding1_2"]
    # Quad1:
    #     Selects first-order quadrupolar frequency contributions. Equivalent to
    #     ["Quad1_2"]
    # Quad2:
    #     Selects second-order quadrupolar frequency contributions. Equivalent to
    #     ["Quad2_0", "Quad2_2", "Quad2_4"]

    # Shielding1: str = "Shielding1"
    Shielding1_0: str = freq_list_all[0]
    Shielding1_2: str = freq_list_all[1]

    # Quad1: str = "Quad1"
    Quad1_2: str = freq_list_all[2]

    # Quad2: str = "Quad2"
    Quad2_0: str = freq_list_all[3]
    Quad2_2: str = freq_list_all[4]
    Quad2_4: str = freq_list_all[5]

    # @validator("Quad1", pre=True, always=True)
    # def validate_Quad1(cls, v, *, values, **kwargs):
    #     values["Quad1_2"] = "Quad1_2"
    #     return None

    # @validator("Quad2", pre=True, always=True)
    # def validate_Quad2(cls, v, *, values, **kwargs):
    #     for item in ["Quad2_0", "Quad2_2", "Quad2_4"]:
    #         values[item] = item
    #     return None

    def json(self) -> dict:
        """Parse the class object to a JSON compliant python dictionary object."""
        return self.value

    def index(self) -> int:
        """Get the index of enumeration relative to freq_list_all."""
        return freq_list_all.index(self.value)


freq_default = [
    "Shielding1_0",
    "Shielding1_2",
    "Quad1_2",
    "Quad2_0",
    "Quad2_2",
    "Quad2_4",
]
default_freq_contrib = [FrequencyEnum(item) for item in freq_default]
