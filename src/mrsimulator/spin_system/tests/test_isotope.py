import pytest
from mrsimulator.spin_system.isotope import get_isotope_dict
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
    assert silicon.efg_to_Cq == 0
    assert silicon.larmor_freq(B0=11.75) == 99.46962016339306
    assert silicon.B0_to_ref_freq(B0=11.75) == 99.3895867929281
    assert silicon.ref_freq_to_B0(ref_freq=99.3895867929281) == 11.75

    proton = Isotope(symbol="1H")
    assert proton.atomic_number == 1
    assert proton.gyromagnetic_ratio == 42.57747920984721
    assert proton.natural_abundance == 99.985
    assert proton.quadrupole_moment == 0.0
    assert proton.spin == 0.5
    assert silicon.efg_to_Cq == 0
    assert proton.larmor_freq(B0=9.40) == -400.2283045725638
    assert proton.B0_to_ref_freq(B0=9.40) == 400.21604182989006
    assert proton.ref_freq_to_B0(ref_freq=400.21604182989006) == 9.40

    nitrogen = Isotope(symbol="14N")
    assert nitrogen.atomic_number == 7
    assert nitrogen.gyromagnetic_ratio == 3.0777058647746447
    assert nitrogen.natural_abundance == 99.634
    assert nitrogen.quadrupole_moment == 0.0193
    assert nitrogen.spin == 1
    assert nitrogen.efg_to_Cq == 4534820.208433115
    assert nitrogen.larmor_freq(B0=18.79) == -57.830093199115574
    assert nitrogen.B0_to_ref_freq(B0=18.79) == 57.81099284148487
    assert nitrogen.ref_freq_to_B0(ref_freq=57.81099284148487) == 18.79

    error = "Isotope symbol `x` not recognized."
    with pytest.raises(Exception, match=error):
        Isotope(symbol="x")

    error = "Isotope symbol `14F` not recognized."
    with pytest.raises(Exception, match=error):
        Isotope(symbol="14F")

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


def test_copy_from():
    error = "`1H` is an immutable registry symbol. Use a different symbol."
    with pytest.raises(Exception, match=error):
        Isotope.register(symbol="1H", gyromagnetic_ratio=20)

    error = "Cannot copy from unrecognized isotope symbol `unknown`"
    with pytest.raises(Exception, match=error):
        Isotope.register(symbol="foo", copy_from="unknown", gyromagnetic_ratio=20)

    Isotope.register(symbol="1H-a", copy_from="1H")
    proton = Isotope(symbol="1H")
    custom_proton = Isotope(symbol="1H-a")

    assert custom_proton.symbol == "1H-a"
    assert custom_proton.spin == proton.spin
    assert custom_proton.spin_multiplicity == proton.spin_multiplicity
    assert custom_proton.gyromagnetic_ratio == proton.gyromagnetic_ratio
    assert custom_proton.quadrupole_moment == proton.quadrupole_moment
    assert custom_proton.natural_abundance == proton.natural_abundance


def test_parse():
    check_iso = Isotope(symbol="custom")

    # Parse an Isotope object, should just return
    iso = Isotope.parse(check_iso)
    assert iso == check_iso

    # Parse in a string
    iso = Isotope.parse("custom")
    assert iso == check_iso

    # Parse in a known dictionary, no overwriting known items
    iso = Isotope.parse({"symbol": "1H"})
    assert iso == Isotope(symbol="1H")

    iso = Isotope.parse(
        {"symbol": "custom", "spin_multiplicity": 4, "gyromagnetic_ratio": -12.3}
    )
    assert iso == check_iso

    # Parse in a dictionary where values for a custom isotope are updated
    iso = Isotope.parse(
        {"symbol": "custom", "spin_multiplicity": 5, "gyromagnetic_ratio": 24.6}
    )
    assert iso.symbol == "custom"
    assert iso.spin == 2  # Value has been updated
    assert iso.spin_multiplicity == 5  # Value has been updated
    assert iso.gyromagnetic_ratio == 24.6  # Value has been updated
    assert iso.quadrupole_moment == 0  # default
    assert iso.natural_abundance == 100  # default

    error = "Cannot parse type <class 'int'> into an Isotope"
    # Too small of integer
    with pytest.raises(Exception, match=error):
        Isotope.parse(123)


def test_validate_isotope_kwargs():
    error = "Isotope spin_multiplicity value must an integer value greater than one."
    # Too small of integer
    with pytest.raises(Exception, match=error):
        Isotope.register(symbol="foo", spin_multiplicity=1)
    # Not an integer
    with pytest.raises(Exception, match=error):
        Isotope.register(symbol="foo", spin_multiplicity=0.2)

    # Natural abundance out of range
    error = "Abundance must be between 0 and 100, inclusive."
    with pytest.raises(Exception, match=error):
        Isotope.register(symbol="foo", natural_abundance=123)
    with pytest.raises(Exception, match=error):
        Isotope.register(symbol="foo", natural_abundance=-10)

    # Atomic number not an integer
    error = "Atomic number must be an integer value."
    with pytest.raises(Exception, match=error):
        Isotope.register(symbol="foo", atomic_number=1.2)


def test_get_isotope_dict():
    error = "Isotope symbol `unknown` not recognized."
    with pytest.raises(Exception, match=error):
        get_isotope_dict(isotope_symbol="unknown")
