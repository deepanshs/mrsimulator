# -*- coding: utf-8 -*-
import sys

import numpy as np


def generate_cos_spherical_distribution(npoints):
    cos_phi = np.cos(2 * np.pi * np.random.rand(npoints))
    cos_theta = 2 * np.random.rand(npoints) - 1
    return cos_phi, cos_theta


def generate_spherical_distribution(npoints):
    phi = 360.0 * np.random.rand(npoints)
    theta = np.arccos(2 * np.random.rand(npoints) - 1) * 180.0 / np.pi
    return phi, theta


def get_CSA_frequencies(three_cos_m1, two_phi2_m1, iso, zeta, eta):
    freq = iso + 0.5 * zeta * (three_cos_m1 + eta * (sin_theta2 * two_phi2_m1))
    return freq


def brute_force_bin(freq, npts, spectral_width):
    delta = spectral_width / npts
    range_ = (
        np.asarray([-spectral_width / 2.0, spectral_width / 2.0]) - delta / 2
    )
    x, y = np.histogram(freq, npts, range_)
    return x, y


def test00(three_cos_m1, two_phi2_m1):
    freq = get_CSA_frequencies(
        three_cos_m1, two_phi2_m1, iso=0, zeta=8000, eta=0.5
    )
    x, y = brute_force_bin(freq, npts=2048, spectral_width=100000)
    np.savetxt("test00.csv", np.asarray([x]).T, fmt="%.4f", delimiter=",")


def test01(three_cos_m1, two_phi2_m1):
    freq = get_CSA_frequencies(
        three_cos_m1, two_phi2_m1, iso=2000, zeta=6000, eta=0.75
    )
    x, y = brute_force_bin(freq, npts=2048, spectral_width=100000)
    np.savetxt("test01.csv", np.asarray([x]).T, fmt="%.4f", delimiter=",")


def test02(three_cos_m1, two_phi2_m1):
    freq = get_CSA_frequencies(
        three_cos_m1, two_phi2_m1, iso=-120, zeta=9340, eta=0.0
    )
    x, y = brute_force_bin(freq, npts=2048, spectral_width=100000)
    np.savetxt("test02.csv", np.asarray([x]).T, fmt="%.4f", delimiter=",")


def test03(three_cos_m1, two_phi2_m1):
    freq = get_CSA_frequencies(
        three_cos_m1, two_phi2_m1, iso=200, zeta=7310, eta=1.0
    )
    x, y = brute_force_bin(freq, npts=2048, spectral_width=100000)
    np.savetxt("test03.csv", np.asarray([x]).T, fmt="%.4f", delimiter=",")


def test04(three_cos_m1, two_phi2_m1):
    freq = get_CSA_frequencies(
        three_cos_m1, two_phi2_m1, iso=1200, zeta=5310, eta=0.31
    )
    x1, y1 = brute_force_bin(freq, npts=2048, spectral_width=100000)
    freq = get_CSA_frequencies(
        three_cos_m1, two_phi2_m1, iso=0, zeta=-3310, eta=0.86
    )
    x2, y2 = brute_force_bin(freq, npts=2048, spectral_width=100000)
    np.savetxt(
        "test04.csv", np.asarray([x1 + x2]).T, fmt="%.4f", delimiter=","
    )


if __name__ == "__main__":
    npoints = 1 * 1000 * 1000000
    cos_phi, cos_theta = generate_cos_spherical_distribution(npoints)
    cos_theta2 = cos_theta ** 2
    cos_theta = None
    del cos_theta
    sin_theta2 = 1.0 - cos_theta2
    three_cos_m1 = 3.0 * cos_theta2 - 1.0
    cos_theta2 = None
    del cos_theta2
    two_phi2_m1 = 2.0 * cos_phi ** 2 - 1.0
    cos_phi = None
    del cos_phi

    test00(three_cos_m1, two_phi2_m1)
    test01(three_cos_m1, two_phi2_m1)
    test02(three_cos_m1, two_phi2_m1)
    test03(three_cos_m1, two_phi2_m1)
    test04(three_cos_m1, two_phi2_m1)

    sys.exit(0)
