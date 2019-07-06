cdef extern from "angular_momentum.h":
    void __wigner_d_matrix(int l, int n, double *angle, double *wigner)

    void __wigner_d_matrix_cosine(int l, int n, double *cos_angle,
                                  double *wigner)

    void __wigner_rotation(int l, int n, double *wigner, double *cos_alpha,
                           double complex *R_in, double complex *R_out)

    void __wigner_dm0_vector(int l, double beta, double *R_out)

cdef extern from "powder_setup.h":
    void __powder_averaging_setup(
        int nt,
        double *cosAlpha,
        double *cosBeta,
        double *amp,
        int space)   # 1 for octant, 2 for hemisphere and 4 for sphere

cdef extern from "interpolation.h":
    void triangle_interpolation(
        double *freq1,
        double *freq2,
        double *freq3,
        double *amp,
        double *spec,
        int *points)

cdef extern from "octahedron.h":
    void octahedronInterpolation(
        double *spec,
        double *freq,
        int nt,
        double *amp,
        int stride,
        int m)

cdef extern from "spinning_sidebands.h":
    void __get_pre_phase_components(
        int number_of_sidebands,
        double spin_frequency,
        double complex *pre_phase)
