# import json
import filecmp
import json
import os

import pytest
from monty.serialization import loadfn
from mrsimulator import Mrsimulator
from mrsimulator import parse
from mrsimulator import update_old_dict_struct
from mrsimulator import update_old_file_struct
from mrsimulator.utils.error import FileConversionError


__author__ = "Matthew D. Giammar"
__email__ = "giammar.7@osu.edu"


MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
test_data = loadfn(os.path.join(MODULE_DIR, "test_data.json"))


def test_parse_old_structure():
    should_be = test_data["serialization_should_be"]
    sim, sp, app = parse(should_be)
    e = ".*An older JSON structure was detected.*"

    # Test old root level struct
    old_struct = test_data["old_root_level"]
    with pytest.warns(Warning, match=e):
        test_sim, test_sp, test_app = parse(old_struct)

    assert test_sim == sim
    assert test_sp == sp
    assert test_app == app

    # Test old root level and convert transition queries
    old_queries = test_data["old_root_level_and_transition_queries"]
    with pytest.warns(Warning, match=e):
        test_sim, test_sp, test_app = parse(old_queries)

    assert test_sim == sim
    assert test_sp == sp
    assert test_app == app

    # Test for error raised when parsing bad structure
    bad_struct = test_data["old_root_level"]
    bad_struct["spin_systems"][0]["bad_key"] = "I will raise an exception"
    with pytest.raises(FileConversionError):
        parse(bad_struct)


def test_update_old_dict_struct():
    old_struct = test_data["old_root_level_and_transition_queries"]
    should_be = test_data["serialization_should_be"]
    assert update_old_dict_struct(old_struct) == Mrsimulator.parse(should_be).json()


def test_update_old_file_struct():
    old_struct = test_data["old_root_level_and_transition_queries"]

    # Save old dictionary to mrsim file
    with open("old_struct.mrsim", "w", encoding="utf8") as outfile:
        json.dump(old_struct, outfile, separators=(",", ":"))

    # Load new structure and save as file to compare later
    should_be = test_data["serialization_should_be"]
    Mrsimulator().parse(should_be).save("should_be.mrsim")

    # Convert and save old_struct.mrsim to converted.mrsim
    old_struct_converted = update_old_file_struct("old_struct.mrsim", "converted.mrsim")

    assert old_struct_converted == Mrsimulator.parse(should_be).json()
    assert filecmp.cmp("converted.mrsim", "should_be.mrsim")

    os.remove("old_struct.mrsim")
    os.remove("should_be.mrsim")
    os.remove("converted.mrsim")
