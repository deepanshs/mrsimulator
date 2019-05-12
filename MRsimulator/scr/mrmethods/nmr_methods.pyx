cimport nmr_methods as clib
cimport numpy as np
import numpy as np
# from nmr.lib import EulerAnglesInRadians
import cython
from .utils import __get_spin_attribute__

@cython.boundscheck(False)
@cython.wraparound(False)
def zg(dict spectrometer,
       list sites,
       dict observed,
       int number_of_sidebands=64,
       int averaging_size=128):
    """
    Parameters
    ----------
    spectrometer:
        A dictionary with spectrometer parameters.
    sites:
        An array of nuclear sites. Each nuclear site is described by 
        a site dictionary. See the description below.
    observed:
        A dictionary of obsevable parameters.

    Description
    -----------
    A spectrometer dictionary consists of the following keywords:

    ``magnetic_flux_density``: A float describing the magnetic flux density
                                 of the spectromenter (in Tesla).
    ``sample_rotation_frequency``: A float describing the sample rotation
                                 frequency (in Hz).
    ``sample_rotation_axis``: A dictionary describing the rotation axis with
                                the following keywords-
        ``polar_angle``: A float describing the angle (in degrees) of the
                           rotating sample with respect to the lab frame z-axis.
        ``azimuthal_angle``: A float describing the angle (in degrees) of the
                           rotating sample with respect to the lab frame x-axis.
    
    An example of the spectrometer dictionary is,

        >>> "spectrometer": {
        ...     "magnetic_flux_density": 9.4,
        ...     "sample_rotation_frequency": 1000,
        ...     "sample_rotation_axis": {
        ...         "polar_angle": 54.735,
        ...         "azimuthal_angle": 0
        ...     }
        ... }
    
    which describes a sample spinning at the magic angle (54.745 degree) at
    1000 Hz inside a 9.4 T NMR spectromenter.

    ---
    A site dictionary consists of the following keywords:

    ``isotope_symbol``: A string with the isotope symbol of the nucleus, eg.
                          '1H' or '29Si'
    ``isotropic_chemical_shift``: A float with the isotropic chemical shift
                          (in Hz).
    ``shielding_symmetric``: A python dictionary describing the spin contributions
                          of a second rank symmetric tensor with the following keywords
        ``anisotropy``: A float describing the strength of the anisotropy (in Hz).
        ``asymmetry``: A float describing the asymmetry of the second rank symmetric
                         tensor. This value ranges from [0,1].

    An example of the site dictionary is
    
        >>> {
        ...     "isotope_symbol": "13C",
        ...     "isotropic_chemical_shift": 100,
        ...     "shielding_symmetric": {
        ...         "anisotropy": 12000,
        ...         "asymmetry": 0.7
        ...     }
        ... }
    
    which describes a '13C' nucleus with an isotropic chemical shift of 100 Hz
    and a second rank symmetric tensor with a 12000 Hz anisotropy and an asymmetry
    of 0.7.

    ---
    The above example 
    """
    # spectrometer
    cdef double sample_rotation_frequency = spectrometer['sample_rotation_frequency']
    cdef double rotor_angle = spectrometer['sample_rotation_axis']['polar_angle']
    B0 = spectrometer['magnetic_flux_density']

    # print ('\nSetting up the virtual NMR spectrometer')
    # print ('---------------------------------------')
    # print (f'Adjusting the virtual magnetic flux density to {B0} T')
    # print (f'Setting sample rotation angle to {rotor_angle} degree')
    # print (f'Setting sample rotation frequency to {sample_rotation_frequency} Hz')
    rotor_angle *= np.pi/180.0

    # spin observed -----------------------------------------------
    obs_spin = observed['isotope_symbol']
    obs_dict = __get_spin_attribute__[obs_spin]

    cdef double spin_quantum_number = obs_dict['spin']
    cdef double larmor_frequency = obs_dict['gyromagnetic_ratio']*B0

    cdef int number_of_points = observed['number_of_points']
    cdef double frequency_bandwidth = observed['frequency_bandwidth']
    cdef double increment = frequency_bandwidth/number_of_points
    cdef double reference_offset = observed['reference_offset'] - frequency_bandwidth/2.0

    # print ((f'Detecting {obs_spin}(I={spin_quantum_number}) isotope '
    #         f'with precession frequency {larmor_frequency} MHz'))
    # print ((f'Recording {obs_spin} spectrum with {number_of_points} '
    #         f'points over a {frequency_bandwidth} Hz bandwidth with a '
    #         f'reference offset of {reference_offset} Hz.'))

    cdef np.ndarray[double, ndim=1] transition_array = \
                        np.asarray(observed['spin_transitions'], dtype=np.float64).ravel()
    cdef int number_of_transitions = transition_array.size/2

    # ---------------------------------------------------------------------
    # sites
    sub_sites = [site for site in sites if site['isotope_symbol'] == obs_spin]

    cdef int number_of_sites= len(sub_sites)

    # site specification
    # CSA
    cdef np.ndarray[double] iso_n = np.empty(number_of_sites, dtype=np.float64)
    cdef np.ndarray[double] aniso_n = np.empty(number_of_sites, dtype=np.float64)
    cdef np.ndarray[double] eta_n = np.empty(number_of_sites, dtype=np.float64)
    
    # quad
    cdef np.ndarray[double] Cq_e = np.zeros(number_of_sites, dtype=np.float64)
    cdef np.ndarray[double] eta_e = np.zeros(number_of_sites, dtype=np.float64)

    cdef np.ndarray[double] D_c = np.zeros(number_of_sites, dtype=np.float64)

    # cdef int i, size = len(sites)
    cdef int i
    for i in range(number_of_sites):
        site = sub_sites[i]

        # CSA tensor
        iso_n[i] = site['isotropic_chemical_shift']
        aniso_n[i] = site['shielding_symmetric']['anisotropy']
        eta_n[i] = site['shielding_symmetric']['asymmetry']
        # print(f'\n{obs_spin} site 1\n----------------------------')
        # print(f'isotropic chemical shift = {iso_n[i]} Hz')
        # print(f'chemical shift anisotropy = {aniso_n[i]} Hz')
        # print(f'chemical shift asymmetry = {eta_n[i]} Hz')
        
        # quad tensor
        # if spin.electric_quadrupole_tensor is not ():
        #     Cq_e[i] = spin.electric_quadrupole_tensor.Cq.to('Hz').value
        #     eta_e[i] = spin.electric_quadrupole_tensor.eta

    # print('\n')
    # print(transition_array)
    # print(number_of_transitions)

    # cdef double m_initial, m_final
    cdef np.ndarray[double, ndim=1] amp = np.zeros(number_of_points * number_of_sites) 
    cdef double cpu_time_

    for i in range(number_of_transitions):
        # spin transitions
        # m_initial = transition_array[i][0]
        # m_final   = transition_array[i][1]
        clib.spinning_sideband_core(
                # spectrum information and related amplitude
                &amp[0],
                &cpu_time_,
                reference_offset,
                increment,
                number_of_points,

                spin_quantum_number,
                larmor_frequency,

                # CSA tensor information
                &iso_n[0],
                &aniso_n[0],
                &eta_n[0],

                # quad tensor information
                &Cq_e[0],
                &eta_e[0],
                1, # quad_second_order_c,

                # Dipolar couplings 
                &D_c[0],

                # spin rate, spin angle and number spinning sidebands
                number_of_sidebands,
                sample_rotation_frequency,
                rotor_angle,

                &transition_array[2*i],
                # Euler angle -> principal to molecular frame
                # &omega_PM_c[0],

                # Euler angles for powder averaging scheme
                averaging_size,
                number_of_sites
                )

    freq = np.arange(number_of_points)*increment + reference_offset   

    return freq, amp.reshape(number_of_sites, number_of_points).sum(axis=0), cpu_time_


