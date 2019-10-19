cimport sandbox as clib

cimport numpy as np
import numpy as np
import cython

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]

# @cython.boundscheck(False)
# @cython.wraparound(False)
# def MRS_plan(int geodesic_polyhedron_frequency,
#         int number_of_sidebands)

cdef class averagingScheme:
    cdef clib.MRS_averaging_scheme *scheme

    def __init__(self, nt, allow_fourth_rank=1):
        self.scheme = clib.MRS_create_averaging_scheme(
                                    nt, allow_fourth_rank)

    @property
    def total_orientations(self):
        return self.scheme.total_orientations

    # def set_scheme(scheme, nt, allow_fourth_rank=1):
    #     if scheme == 'octahedron':
    #         self.scheme = clib.MRS_create_averaging_scheme(
    #                                 nt, allow_fourth_rank)

# cdef class MRSPlan:
#     cdef clib.MRS_plan *plan

#     def __init__(self, averaging_scheme, number_of_sidebands,
#                 sample_rotation_frequency_in_Hz, rotor_angle_in_rad,
#                 increment, allow_fourth_rank):
#         self.averaging_scheme = averaging_scheme
#         self.plan = clib.MRS_create_plan(averaging_scheme.scheme,
#                         number_of_sidebands, sample_rotation_frequency_in_Hz,
#                         rotor_angle_in_rad, increment, allow_fourth_rank)

#     @property
#     def number_of_sidebands(self):
#         return self.plan.number_of_sidebands

#     @property
#     def sample_rotation_frequency_in_Hz(self):
#         return self.plan.sample_rotation_frequency_in_Hz

#     @property
#     def rotor_angle_in_rad(self):
#         return self.plan.rotor_angle_in_rad

#     @property
#     def increment(self):
#         return self.plan.increment

#     @property
#     def allow_fourth_rank(self):
#         return self.plan.allow_fourth_rank

#     def evaluate(self, R0, R2, R4):
#         cdef np.ndarray[double complex] R2_c = np.asarray(R2, dtype=np.complex128)
#         cdef np.ndarray[double complex] R4_c = np.asarray(R4, dtype=np.complex128)
#         clib.MRS_get_amplitudes_from_plan(self.plan, &R2_c[0], &R4_c[0])
#         clib.MRS_get_frequencies_from_plan(self.plan, R0)
#         # side_band = np.
#         return self.plan.vector[::2].reshape(self.plan.number_of_sidebands,
#                                     self.plan.total_orientations)
