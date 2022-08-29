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

    with pytest.raises(Exception, match="Could not parse isotope string"):
        Isotope(symbol="x")

    with pytest.raises(Exception, match="Could not parse isotope string"):
        Isotope(symbol="14F")

    assert nitrogen.json() == "14N"


def test_custom_isotope():
    custom_isotope_dict = {
        "symbol": "custom isotope",
        "spin": 4,
        "natural_abundance": 45.6,
        "gyromagnetic_ratio": -8,
        "quadrupole_moment": 0.123,
        "atomic_number": 0,
    }

    iso = Isotope.add_new(**custom_isotope_dict)

    assert iso.atomic_number == 0
    assert iso.gyromagnetic_ratio == -8
    assert iso.natural_abundance == 45.6
    assert iso.quadrupole_moment == 0.123
    assert iso.spin == 4
    assert iso.json() == custom_isotope_dict

    # Creating new isotope with same symbol
    new_iso = Isotope.add_new(symbol="custom isotope", spin=4, gyromagnetic_ratio=-8)

    assert new_iso.json() == iso.json()

    # Error on creating new isotope with overlapping symbol
    with pytest.raises(ValueError, match=".*cannot match a real isotope symbol.*"):
        Isotope.add_new(symbol="1H", spin=0.5, gyromagnetic_ratio=42.57747920984721)

    # Error on spin less than zero or not half-integer
    with pytest.raises(ValueError, match=".*Isotope spin value must be greater than.*"):
        Isotope.add_new(symbol="bad iso", spin=0, gyromagnetic_ratio=10)

    with pytest.raises(ValueError, match=".*Isotope spin value must be greater than.*"):
        Isotope.add_new(symbol="bad iso", spin=1.75, gyromagnetic_ratio=10)

    # Error on abundance outsize of [0, 100]
    with pytest.raises(ValueError, match=".*Abundance must be between 0 and 100.*"):
        Isotope.add_new(
            symbol="bad iso",
            spin=1.5,
            gyromagnetic_ratio=10,
            natural_abundance=150,
        )
