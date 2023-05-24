import pytest
from mrsimulator.spin_system.isotope import Isotope

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


def test_isotope():
    silicon = Isotope(symbol="29Si")

    assert silicon.symbol == "29Si"
    assert silicon.atomic_number == 14
    assert silicon.gyromagnetic_ratio == -8.465499588373877
    assert silicon.natural_abundance == 4.683
    assert silicon.quadrupole_moment == 0.0
    assert silicon.spin == 0.5
    assert silicon.larmor_freq(B0=11.75) == 99.46962016339306

    proton = Isotope(symbol="1H")
    assert proton.atomic_number == 1
    assert proton.gyromagnetic_ratio == 42.57747920984721
    assert proton.natural_abundance == 99.985
    assert proton.quadrupole_moment == 0.0
    assert proton.spin == 0.5
    assert proton.larmor_freq(B0=9.40) == -400.2283045725638

    nitrogen = Isotope(symbol="14N")
    assert nitrogen.atomic_number == 7
    assert nitrogen.gyromagnetic_ratio == 3.0777058647746447
    assert nitrogen.natural_abundance == 99.634
    assert nitrogen.quadrupole_moment == 0.0193
    assert nitrogen.spin == 1
    assert nitrogen.larmor_freq(B0=18.79) == -57.830093199115574

    error = "is immutable and does not support item assignment"
    with pytest.raises(Exception, match=error):
        nitrogen.spin_multiplicity = 4

    error = "spin_multiplicity for 1H cannot be assigned."
    with pytest.raises(Exception, match=error):
        _ = Isotope(symbol="1H", spin_multiplicity=5)

    assert nitrogen.json() == "14N"


def test_custom_isotope():
    Isotope.register(
        symbol="custom",
        spin_multiplicity=4,
        gyromagnetic_ratio=-12.3,
        quadrupole_moment=0.1,
        natural_abundance=50,
    )

    custom = Isotope(symbol="custom")

    assert isinstance(custom, Isotope)
    assert custom.symbol == "custom"
    assert custom.spin == 1.5
    assert custom.spin_multiplicity == 4
    assert custom.gyromagnetic_ratio == -12.3
    assert custom.quadrupole_moment == 0.1
    assert custom.natural_abundance == 50

    assert custom.json() == {
        "symbol": "custom",
        "spin_multiplicity": 4,
        "gyromagnetic_ratio": -12.3,
        "quadrupole_moment": 0.1,
        "natural_abundance": 50,
        "atomic_number": 0,
    }
