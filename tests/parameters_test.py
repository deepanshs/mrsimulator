# -*- coding: utf-8 -*-
"""Parameter test"""
import csdmpy as cp
import numpy as np
from lmfit import Parameters
from mrsimulator import Dimension
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SymmetricTensor as st
from mrsimulator.spectral_fitting import make_fitting_parameters


sim = Simulator()

H = Site(
    isotope="1H", isotropic_chemical_shift=10, shielding_symmetric=st(zeta=5, eta=0.1)
)
isotopomer = Isotopomer(name="H1", sites=[H], abundance=100)

dimension = Dimension(
    isotope="1H",
    magnetic_flux_density=9.4,  # in T
    number_of_points=2046,
    spectral_width=20000,  # in Hz
    reference_offset=0,
    rotor_frequency=14000,  # in Hz
)

sim.isotopomers += [isotopomer]
sim.dimensions += [dimension]

params = make_fitting_parameters(sim)


test_params = Parameters()
test_params.add(name="isotopomers9109346sites9109346isotropic_chemical_shift", value=10)
test_params.add(name="isotopomers9109346sites9109346shielding_symmetric46zeta", value=5)
test_params.add(
    name="isotopomers9109346sites9109346shielding_symmetric46eta",
    value=0.1,
    min=0,
    max=1,
)
test_params.add(
    name="isotopomers9109346abundance",
    value=100,
    min=0,
    max=100,
    vary=False,
    expr="100",
)


def param_test():
    assert params.valuesdict() == test_params.valuesdict(), "Parameter creation failed"
