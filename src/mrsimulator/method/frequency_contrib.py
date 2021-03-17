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
        Select first-order and zeroth-rank nuclear shielding frequency contributions.

    Shielding1_2:
        Select first-order and second-rank nuclear shielding frequency contributions.

    Quad1_2:
        Select first-order and second-rank quadrupolar frequency contributions.

    Quad2_0:
        Select second-order and zeroth-rank quadrupolar frequency contributions.

    Quad2_2:
        Select second-order and second-rank quadrupolar frequency contributions.

    Quad2_4:
        Select second-order and fourth-rank quadrupolar frequency contributions.
    """

    Shielding1_0 = freq_list_all[0]
    Shielding1_2 = freq_list_all[1]
    Quad1_2 = freq_list_all[2]
    Quad2_0 = freq_list_all[3]
    Quad2_2 = freq_list_all[4]
    Quad2_4 = freq_list_all[5]

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
