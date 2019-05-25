

def csa_static():
    temp = {
        "isotopomers": [
            {
                "sites": [
                    {
                        "isotope_symbol": "1H",
                        "isotropic_chemical_shift": "0 Hz",
                        "shielding_symmetric": {
                            "anisotropy": "-3.89 kHz",
                            "asymmetry": 0.25
                        }
                    }
                ],
                "abundance": "100 %"
            }
        ],
        "spectrum": {
            "direct_dimension":{
                "magnetic_flux_density": "9.4 T",
                "rotor_frequency": "0 kHz",
                "rotor_angle": "54.735 deg",
                "number_of_points": 2048,
                "spectral_width": "25 kHz",
                "reference_offset": "0 Hz",
                "nucleus": "1H"
            }
        }
    }
    return temp["isotopomers"], temp["spectrum"]

def csa_mas():
    temp = {
        "isotopomers": [
            {
                "sites": [
                    {
                        "isotope_symbol": "13C",
                        "isotropic_chemical_shift": "0 Hz",
                        "shielding_symmetric": {
                            "anisotropy": "12.89 kHz",
                            "asymmetry": 0.5
                        }
                    }
                ],
                "abundance": "100 %"
            }
        ],
        "spectrum": {
            "direct_dimension":{
                "magnetic_flux_density": "9.4 T",
                "rotor_frequency": "5 kHz",
                "rotor_angle": "54.735 deg",
                "number_of_points": 2048,
                "spectral_width": "100 kHz",
                "reference_offset": "0 Hz",
                "nucleus": "13C"
            }
        }
    }
    return temp["isotopomers"], temp["spectrum"]