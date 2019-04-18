
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

cdef extern from "spinning_sidebands.h":
    void lineshape_cas_spinning_sideband_core(
        # spectrum information and related amplitude
        double * spec,
        double * cpu_time_,
        double spectral_start,
        double spectral_increment,
        int number_of_points,

        # Pointer to spin quantum numbers for every spin.
        double *qunatum_number,
        double *wo,

        # Pointer to the array of CSA tensor information in the PAS for every spin. 
        double *iso,
        double *aniso,
        double *eta,

        # Pointer to the array of quadrupole tensor information in the PAS for every spin. 
        double *Cq_e,                       # The Cq of the quadrupole center.
        double *eta_e,                      # The asymmetry term of the tensor.
        int quadSecondOrder,                # Quad theory for second order, 

        # spin rate, spin angle and number spinning sidebands
        int ph_step,
        double spin_frequency,
        double rotor_angle,

        # The transition as transition[0] = mi and transition[1] = mf
        double *transition,     
                
        # Euler angle -> principal to molecular frame
        # double *omega_PM,

        # Euler angles for powder averaging scheme
        int nt,
        unsigned int number_of_site)
    #     int averaging_scheme,
        
    # )

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

# cdef extern from "MRAngularMomentum.h":
#     void fullWigner_d(double *wigner, double l, double beta)
#     double wigner_d(double l, double m1, double m2, double beta)
