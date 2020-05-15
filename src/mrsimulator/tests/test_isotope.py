# -*- coding: utf-8 -*-
import pytest
from mrsimulator.isotope import Isotope


def test_isotope():
    silicon = Isotope(symbol="29Si")

    assert silicon.symbol == "29Si"
    assert silicon.atomic_number == 14
    assert silicon.gyromagnetic_ratio == -8.465499
    assert silicon.natural_abundance == 4.683
    assert silicon.quadrupole_moment == 0.0
    assert silicon.spin == 0.5

    proton = Isotope(symbol="1H")
    assert proton.atomic_number == 1
    assert proton.gyromagnetic_ratio == 42.57748
    assert proton.natural_abundance == 99.985
    assert proton.quadrupole_moment == 0.0
    assert proton.spin == 0.5

    nitrogen = Isotope(symbol="14N")
    assert nitrogen.atomic_number == 7
    assert nitrogen.gyromagnetic_ratio == 3.077706
    assert nitrogen.natural_abundance == 99.634
    assert nitrogen.quadrupole_moment == 0.0193
    assert nitrogen.spin == 1

    with pytest.raises(Exception, match="Could not parse isotope string"):
        Isotope(symbol="x")

    with pytest.raises(Exception, match="Could not parse isotope string"):
        Isotope(symbol="14F")

    assert nitrogen.to_dict_with_units() == "14N"
