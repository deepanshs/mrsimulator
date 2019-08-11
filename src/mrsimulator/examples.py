# -*- coding: utf-8 -*-
__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


def csa_static():
    temp = {
        "isotopomers": [
            {
                "sites": [
                    {
                        "isotope": "1H",
                        "isotropic_chemical_shift": "0 ppm",
                        "shielding_symmetric": {
                            "anisotropy": "13.89 ppm",
                            "asymmetry": 0.25,
                        },
                    }
                ],
                "abundance": "100 %",
            }
        ],
        "spectrum": {
            "direct_dimension": {
                "magnetic_flux_density": "9.4 T",
                "rotor_frequency": "0 kHz",
                "rotor_angle": "54.735 deg",
                "number_of_points": 2048,
                "spectral_width": "25 kHz",
                "reference_offset": "0 Hz",
                "isotope": "1H",
            }
        },
    }
    return temp["isotopomers"], temp["spectrum"]


def csa_mas():
    temp = {
        "isotopomers": [
            {
                "sites": [
                    {
                        "isotope": "1H",
                        "isotropic_chemical_shift": "0 ppm",
                        "shielding_symmetric": {
                            "anisotropy": "13.89 ppm",
                            "asymmetry": 0.25,
                        },
                    }
                ],
                "abundance": "100 %",
            }
        ],
        "spectrum": {
            "direct_dimension": {
                "magnetic_flux_density": "9.4 T",
                "rotor_frequency": "1 kHz",
                "rotor_angle": "54.735 deg",
                "number_of_points": 2048,
                "spectral_width": "25 kHz",
                "reference_offset": "0 Hz",
                "isotope": "1H",
            }
        },
    }
    return temp["isotopomers"], temp["spectrum"]
