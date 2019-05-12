

cdef extern from "csa_static_lineshape.h":
    void lineshape_csa_static(
        double * spec,
        double * cpu_time_,
        int m,
        int nt,
        double fstart,
        double iso,
        double fwidth,
        double aniso,
        double eta,
        int octa,
        int npros
    ) 

cdef extern from "powder_setup.h":
    void rasterization(double * grid,
                   double *v0,
                   double *v1,
                   double *v2,
                   int rows,
                   int columns)

cdef extern from "spinning_sidebands.h":
    void lineshape_cas_spinning_sideband_cython_angles(
        # spectrum information and related amplitude
        double * spec,
        double * cpu_time_,
        double frequency_start,
        double frequency_bandwidth,
        int number_of_points,

        # CSA tensor information
        double iso,
        double aniso,
        double eta,

        # spin rate, spin angle and number spinning sidebands
        int ph_step,
        double spin_frequency,
        double rotor_angle,

        # Euler angle -> principal to molecular frame
        double *omega_PM,

        # Euler angles for powder averaging scheme
        int averaging_scheme,
        int averaging_size)