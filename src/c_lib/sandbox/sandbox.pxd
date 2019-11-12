from libcpp cimport bool as bool_t

cdef extern from "averaging_scheme.h":
    ctypedef struct MRS_averaging_scheme:
        unsigned int total_orientations

    MRS_averaging_scheme * MRS_create_averaging_scheme(
                            unsigned int geodesic_polyhedron_frequency,
                            bool_t allow_fourth_rank,
                            unsigned int integration_volume)

    void MRS_free_averaging_scheme(MRS_averaging_scheme *scheme)

cdef extern from "mrsimulator.h":

    ctypedef struct MRS_plan:
        MRS_averaging_scheme *averaging_scheme
        int number_of_sidebands
        double sample_rotation_frequency_in_Hz
        double rotor_angle_in_rad
        double complex *vector

    MRS_plan *MRS_create_plan(MRS_averaging_scheme *scheme, int number_of_sidebands,
                          double sample_rotation_frequency_in_Hz,
                          double rotor_angle_in_rad, double increment,
                          bool_t allow_fourth_rank)
    void MRS_free_plan(MRS_plan *plan)
    void MRS_get_amplitudes_from_plan(MRS_plan *plan, double complex *R2,
                                  double complex *R4)
    void MRS_get_frequencies_from_plan(MRS_plan *plan, double R0)
