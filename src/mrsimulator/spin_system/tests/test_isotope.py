import pytest
from mrsimulator.spin_system.isotope import get_isotope_data
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
    Isotope.new(symbol="Si_a", copy_from="29Si")
    Isotope.new(symbol="Si_b", copy_from="29Si")

    assert get_isotope_data("Si_a") == get_isotope_data("29Si")
    assert get_isotope_data("Si_b") == get_isotope_data("Si_a")

    # Overlapping isotope symbols
    error = r".*The new symbol for the copied isotope cannot match any previous.*"
    with pytest.raises(ValueError, match=error):
        Isotope.new(symbol="17O", copy_from="17O")

    error = r".*The copy_from symbol must match a known isotope.*"
    with pytest.raises(ValueError, match=error):
        Isotope.new(symbol="new", copy_from="unknown")
