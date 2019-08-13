# -*- coding: utf-8 -*-
import os.path
import pytest
from monty.serialization import loadfn

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_FOLDER = os.path.abspath(os.path.join(MODULE_DIR, "..", "..", "..", "test_files"))


@pytest.fixture
def mas_data():
    return loadfn(os.path.join(TEST_FOLDER, "mas.json"))


@pytest.fixture
def static_data():
    return loadfn(os.path.join(TEST_FOLDER, "static.json"))
