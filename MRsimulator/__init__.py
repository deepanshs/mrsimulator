import numpy as np
from collections import namedtuple
from dataclasses import dataclass
from .unit import string_to_quantity
from astropy.units import Quantity
from astropy import units as u

from . import core, parameterized_tensor


__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.90@osu.edu"


# @dataclass
# class Hamiltonian:
#     _first_order_nuclear_shielding_: int = 0
#     _first_order_electric_quadrupole_coupling_: int =0
#     _second_order_electric_quadrupole_coupling_: int=0

#     def __add__(self, other):
#         i_ = self._first_order_nuclear_shielding_
#         i_+= other._first_order_nuclear_shielding_
#         j_ = self._first_order_electric_quadrupole_coupling_
#         j_+= other._first_order_electric_quadrupole_coupling_
#         k_ = self._second_order_electric_quadrupole_coupling_
#         k_+= other._second_order_electric_quadrupole_coupling_
#         return Hamiltonian(i_, j_, k_)

# class FirstOrderNuclearShielding:
#     def __new__(self):
#         return Hamiltonian(_first_order_nuclear_shielding_ = 1)

# class FirstOrderElectricQuadrupoleCoupling:
#     def __new__(self):
#         return Hamiltonian(_first_order_electric_quadrupole_coupling_ = 1)

# class SecondOrderElectricQuadrupoleCoupling:
#     def __new__(self):
#         return Hamiltonian(_second_order_electric_quadrupole_coupling_ = 1)













