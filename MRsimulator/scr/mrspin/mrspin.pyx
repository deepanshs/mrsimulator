
cimport numpy as np

cdef struct c_NuclearShieldingTensor:
    double isotropicNuclearShiftInPPM          # Isotropic chemical shift (ppm)                   #
    double shieldingTensorAnisotropyInPPM      # Nuclear shielding anisotropy (ppm)               #
    double shieldingTensorAsymmetry            # Nuclear shielding asymmetry parameter            #
    double *shieldingTensorOrientation         # Nuclear shielding PAS to CRS euler angles (rad.) #

cdef struct QuadrupolarTensor:
    double quadrupolarTensorConstantInHz       # Quadrupolar coupling constant (Hz)               #
    double quadrupolarTensorAsymmetry          # Quadrupolar asymmetry parameter                  #
    double *quadrupolarTensorOrientation       # Quadrupolar PAS to CRS euler angles (rad.)       #

cdef struct c_MRsite:
    float spinQuantumNumber
    NuclearShieldingTensor nuclearShieldingInPAS
    NuclearShieldingTensor quandupolarCouplingInPAS

cdef struct transition:
    np.ndarray m_initial
    double m_final




