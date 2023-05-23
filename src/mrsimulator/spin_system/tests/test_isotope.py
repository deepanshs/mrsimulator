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

    # Full isotope dictionary serialized
    assert nitrogen.json() == {
        "spin_multiplicity": 3,
        "natural_abundance": 99.634,
        "gyromagnetic_ratio": 3.0777058647746447,
        "quadrupole_moment": 0.0193,
        "atomic_number": 7,
        "symbol": "14N",
    }


def test_custom_isotope():
    custom = Isotope(
        symbol="custom",
        spin_multiplicity=4,
        gyromagnetic_ratio=-12.3,
        quadrupole_moment=0.1,
        natural_abundance=50,
    )

    assert isinstance(custom, Isotope)
    assert custom.symbol == "custom"
    assert custom.spin == 1.5
    assert custom.spin_multiplicity == 4
    assert custom.gyromagnetic_ratio == -12.3
    assert custom.quadrupole_moment == 0.1
    assert custom.natural_abundance == 50

    # Instantiating Isotope from just string
    custom_from_string = Isotope(symbol="custom")

    assert custom == custom_from_string

    assert custom.json() == {
        "symbol": "custom",
        "spin_multiplicity": 4,
        "gyromagnetic_ratio": -12.3,
        "quadrupole_moment": 0.1,
        "natural_abundance": 50,
        "atomic_number": 0,
    }

    # Extra fields passed to Isotope
    e = "Extra keywords, {'im_an_extra_key'}, are prohibited"
    with pytest.raises(Exception, match=e):
        Isotope(symbol="foo", im_an_extra_key=1)


def test_validate():
    # # Ensuring no overlap on symbols
    # e = "already attributed to another Isotope. All Isotope symbols must be unique;"
    # with pytest.raises(Exception, match=e):
    #     Isotope(symbol="custom", spin_multiplicity=3, gyromagnetic_ratio=1)

    # Errors for arguments out of allowed ranges:
    e = "Isotope spin_multiplicity value must an integer value greater than one"
    with pytest.raises(Exception, match=e):
        Isotope(symbol="foo", spin_multiplicity=1, gyromagnetic_ratio=1)
    with pytest.raises(Exception, match=e):
        Isotope(symbol="foo", spin_multiplicity=1.2, gyromagnetic_ratio=1)

    e = "Abundance must be between 0 and 100, inclusive."
    with pytest.raises(Exception, match=e):
        Isotope(
            symbol="foo",
            spin_multiplicity=3,
            gyromagnetic_ratio=1,
            natural_abundance=150,
        )

    e = "Atomic number must be an integer value."
    with pytest.raises(Exception, match=e):
        Isotope(
            symbol="foo", spin_multiplicity=3, gyromagnetic_ratio=1, atomic_number=5.5
        )


def test_get_isotope():
    iso_1H = Isotope(symbol="1H")
    good_1H_dict = {
        "natural_abundance": 99.985,
        "gyromagnetic_ratio": 42.57747920984721,
        "quadrupole_moment": 0.0,
        "atomic_number": 1,
        "spin_multiplicity": 2,
        "symbol": "1H",
    }
    bad_1H_dict = {
        "natural_abundance": 99.985,
        "gyromagnetic_ratio": 42.57747920984721,
        "quadrupole_moment": 0.0,
        "atomic_number": 1,
        "spin_multiplicity": 3,  # Changed spin multiplicity number
        "symbol": "1H",
    }

    assert Isotope.get_isotope(iso_1H) == iso_1H
    assert Isotope.get_isotope("1H") == iso_1H
    assert Isotope.get_isotope(good_1H_dict) == iso_1H

    # Passing dictionary with different stored data
    e = "cannot be assigned"
    with pytest.raises(Exception, match=e):
        Isotope.get_isotope(val=bad_1H_dict)

    # Passing an unknown symbol to get_isotope
    e = "Could not parse isotope string"
    with pytest.raises(Exception, match=e):
        Isotope.get_isotope(val="unknown_symbol")

    new_custom_isotope = {
        "natural_abundance": 12.34,
        "gyromagnetic_ratio": 5.67,
        "quadrupole_moment": 8.9,
        "atomic_number": -1,
        "spin_multiplicity": 5,
        "symbol": "new_custom_isotope",
    }

    new_iso = Isotope.get_isotope(new_custom_isotope)

    assert new_iso.symbol == "new_custom_isotope"
    assert new_iso.spin == 2
    assert new_iso.gyromagnetic_ratio == 5.67
    assert new_iso.quadrupole_moment == 8.9
    assert new_iso.natural_abundance == 12.34

    # Passing an object other than string, dictionary, or Isotope
    e = "is invalid for this method"
    with pytest.raises(Exception, match=e):
        Isotope.get_isotope(val=10)
