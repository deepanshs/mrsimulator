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
    Shielding1_0 = freq_list_all[0]
    Shielding1_2 = freq_list_all[1]
    Quad1_2 = freq_list_all[2]
    Quad2_0 = freq_list_all[3]
    Quad2_2 = freq_list_all[4]
    Quad2_4 = freq_list_all[5]

    def to_dict_with_units(self):
        return self.value


freq_default = [
    "Shielding1_0",
    "Shielding1_2",
    "Quad1_2",
    "Quad2_0",
    "Quad2_2",
    "Quad2_4",
]
default_freq_contrib = [FrequencyEnum(item) for item in freq_default]
