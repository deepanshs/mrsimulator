cdef extern from "angular_momentum.h":
    void __wigner_d_matrix_cosine(int l, int n, double *cos_angle,
                                  double *wigner)

cdef extern from "powder_setup.h":
    void __powder_averaging_setup(int nt, double *cosAlpha, double *cosBeta,
                                  double *amp, int space)

# cdef extern from "site.h":
#     ctypedef struct site:
#         unsigned int index
#         double isotropic_chemical_shift_in_Hz # Isotropic chemical shift (Hz)
#         double shielding_symmetric_anisotropy_in_Hz     # Nuclear shielding anisotropy (Hz)
#         double shielding_symmetric_asymmetry            # Nuclear shielding asymmetry
#         double shielding_symmetric_orientation[3]          # Nuclear shielding PAS to CRS euler angles (rad.)
#         double quadrupolar_constant_in_Hz     # Quadrupolar coupling constant (Hz)
#         double quadrupolar_asymmetry          # Quadrupolar asymmetry parameter
#         double quadrupolar_orientation[3]        # Quadrupolar PAS to CRS euler angles (rad.)
#         double dipolar_coupling              # dipolar coupling sof the site


cdef extern from "isotopomer_ravel.h":
    ctypedef struct isotopomer_ravel:
        int number_of_sites                    # Number of sites
        float spin                             # The spin quantum number
        double larmor_frequency                # Larmor frequency (MHz)
        double *isotropic_chemical_shift_in_Hz # Isotropic chemical shift (Hz)
        double *shielding_anisotropy_in_Hz     # Nuclear shielding anisotropy (Hz)
        double *shielding_asymmetry            # Nuclear shielding asymmetry
        double *shielding_orientation          # Nuclear shielding PAS to CRS euler angles (rad.)
        double *quadrupolar_constant_in_Hz     # Quadrupolar coupling constant (Hz)
        double *quadrupolar_asymmetry          # Quadrupolar asymmetry parameter
        double *quadrupolar_orientation        # Quadrupolar PAS to CRS euler angles (rad.)
        double *dipolar_couplings              # dipolar coupling stored as list of lists

cdef extern from "spinning_sidebands.h":
    void spinning_sideband_core(
        # spectrum information and related amplitude
        double * spec,
        double spectral_start,
        double spectral_increment,
        int number_of_points,

        isotopomer_ravel *ravel_isotopomer,

        int quad_second_order,                    # Quad theory for second order,
        int remove_second_order_quad_isotropic,   # remove the isotropic contribution from the
                                                  # second order quad Hamiltonian.

        # spin rate, spin angle and number spinning sidebands
        int number_of_sidebands,
        double sample_rotation_frequency_in_Hz,
        double rotor_angle_in_rad,

        # The transition as transition[0] = mi and transition[1] = mf
        double *transition,
        int geodesic_polyhedron_frequency)
