

import numpy as np

class Hamiltonian:
    def __init__(self):
        self.__first_order_nuclear_shielding__ = 0
        self.__first_order_electric_quadrupole_coupling__ = 0
        self.__second_order_electric_quadrupole_coupling__ = 0

class FirstOrderNuclearShielding(Hamiltonian):
    def __init__(self):
        self.__first_order_nuclear_shielding__ = 1

class FirstOrderElectricQuadrupoleCoupling(Hamiltonian):
    def __init__(self):
        self.__first_order_electric_quadrupole_coupling__ = 1

class SecondOrderElectricQuadrupoleCoupling(Hamiltonian):
    def __init__(self):
        self.__second_order_electric_quadrupole_coupling__ = 1

