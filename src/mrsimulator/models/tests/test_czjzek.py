# -*- coding: utf-8 -*-
from os import path

import numpy as np
from mrsimulator.models import extended_czjzek_distribution

# from mrsimulator.models import czjzek_distribution

MODULE_DIR = path.dirname(path.abspath(__file__))


def test_extended_czjzek_eta_distribution_1():
    filename = path.join(MODULE_DIR, "test_data", "eps=0.05.npy")
    with open(filename, "rb") as f:
        data = np.load(f)

    print(data.shape)
    _, eta1 = extended_czjzek_distribution(1, 0.1, eps=0.05, n=np.int(1e6))
    hist1, _ = np.histogram(eta1, bins=100, range=[0, 1])

    message = "failed to compare values with file eps=0.05.npy"
    np.testing.assert_almost_equal(hist1 / 1e6, data[0], decimal=3, err_msg=message)


def test_extended_czjzek_eta_distribution_2():
    filename = path.join(MODULE_DIR, "test_data", "eps=0.2.npy")
    with open(filename, "rb") as f:
        data = np.load(f)

    print(data.shape)
    _, eta1 = extended_czjzek_distribution(1, 0.3, eps=0.2, n=np.int(2.5e6))
    hist1, _ = np.histogram(eta1, bins=100, range=[0, 1])

    message = "failed to compare values with file eps=0.2.npy"
    np.testing.assert_almost_equal(hist1 / 2.5e6, data[1], decimal=3, err_msg=message)


def test_extended_czjzek_eta_distribution_3():
    filename = path.join(MODULE_DIR, "test_data", "eps=0.65.npy")
    with open(filename, "rb") as f:
        data = np.load(f)

    print(data.shape)
    _, eta1 = extended_czjzek_distribution(1, 0.7, eps=0.65, n=np.int(2.5e6))
    hist1, _ = np.histogram(eta1, bins=100, range=[0, 1])

    message = "failed to compare values with file eps=0.05.npy"
    np.testing.assert_almost_equal(hist1 / 2.5e6, data[3], decimal=3, err_msg=message)
