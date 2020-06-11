# -*- coding: utf-8 -*-
"""Parameter test"""
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator.spectral_fitting import make_fitting_parameters
from mrsimulator.tensors import SymmetricTensor as st


sim = Simulator()

H = Site(
    isotope="1H", isotropic_chemical_shift=10, shielding_symmetric=st(zeta=5, eta=0.1)
)
spin_system = SpinSystem(name="H1", sites=[H], abundance=100)

method = BlochDecaySpectrum(
    channels=["1H"],
    magnetic_flux_density=9.4,
    rotor_frequency=14000,
    spectral_dimensions=[
        {"count": 2048, "spectral_width": 20000, "reference_offset": 0}
    ],
)


sim.spin_systems += [spin_system]
sim.methods += [method]

params = make_fitting_parameters(sim)

valuesdict = {
    "ISO_0_SITES_0_isotropic_chemical_shift": 10,
    "ISO_0_SITES_0_shielding_symmetric_zeta": 5,
    "ISO_0_SITES_0_shielding_symmetric_eta": 0.1,
    "ISO_0_abundance": 100,
}


def test_param():
    assert params.valuesdict() == valuesdict, "Parameter creation failed"
