
cdef extern from "spinning_sidebands.h":
    void spinning_sideband_core(
        # spectrum information and related amplitude
        double * spec,
        double * cpu_time_,
        double spectral_start,
        double spectral_increment,
        int number_of_points,

        # Pointer to spin quantum numbers for every spin.
        double spin_quantum_number,
        double larmor_frequency,

        # Pointer to the array of CSA tensor information in the PAS for every spin. 
        double *iso,
        double *aniso,
        double *eta,

        # Pointer to the array of quadrupole tensor information in the PAS for every spin. 
        double *Cq_e,                       # The Cq of the quadrupole center.
        double *eta_e,                      # The asymmetry term of the tensor.
        int quadSecondOrder,                # Quad theory for second order, 

				# Pointer to the array of dipolar tensor information in the PAS. 
        double *D,                          # The dipolar coupling constant.

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