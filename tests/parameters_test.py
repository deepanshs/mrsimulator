# -*- coding: utf-8 -*-
"""Parameter test"""
from lmfit import Parameters
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SymmetricTensor as st
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator.spectral_fitting import make_fitting_parameters


sim = Simulator()

H = Site(
    isotope="1H", isotropic_chemical_shift=10, shielding_symmetric=st(zeta=5, eta=0.1)
)
isotopomer = Isotopomer(name="H1", sites=[H], abundance=100)

method = BlochDecaySpectrum(
    channels=["1H"],
    magnetic_flux_density=9.4,
    rotor_frequency=14000,
    dimensions=[{"count": 2046, "spectral_width": 20000, "reference_offset": 0}],
)


sim.isotopomers += [isotopomer]
sim.methods += [method]

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
