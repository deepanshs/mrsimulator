cimport sandbox as clib
from libcpp cimport bool as bool_t

cimport numpy as np
import numpy as np
import cython

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

__integration_volume_enum__ = {"octant": 0, "hemisphere": 1}
__integration_volume_enum_rev__ = {0: "octant", 1: "hemisphere"}


# cdef class UserDefinedAveragingScheme:
#     cdef clib.MRS_averaging_scheme *scheme

#     def __init__(
#             self,
#             np.ndarray[ndim=1, double] alpha,
#             np.ndarray[ndim=1, double] beta,
#             np.ndarray[ndim=1, double] weight,
#             bool_t allow_fourth_rank=False
#         ):
#         """Create the generic averaging scheme for the simulation.""""
#         if alpha.size != beta.size != weight.size:
#             raise ValueError(
#                 'The length of alpha, beta, and weight array must be equal.'
#             )
#         cdef int n_angles = alpha.size
#         self.scheme = clib.MRS_create_averaging_scheme_from_alpha_beta(&alpha[0],
#                                             &beta[0], &weight[0], n_angles,
#                                             allow_fourth_rank_)
#     @property
#     def interpolation(self):
#         return False

#     @property
#     def total_orientations(self):
#         return self.scheme.total_orientations

#     def __reduce_cython__(self):
#         return 'AveragingScheme'

#     def __eq__(self, other):
#         if isinstance(other, AveragingScheme):
#             check = [
#                 self.alpha == other.alpha,
#                 self.beta == other.beta
#                 self.weight == other.weight
#             ]
#             if np.all(check):
#                 return True
#         return False

#     def __str__(self):
#         return self.__repr__()

#     def __repr__(self):
#         return (
#             "AveragingScheme(total_orientations={0}, interpolation={1})"
#         ).format(self.total_orientations, self.interpolation)

#     # def __del__(self):
#     #     clib.MRS_free_averaging_scheme(self.scheme)


cdef class AveragingScheme:
    cdef clib.MRS_averaging_scheme *scheme
    cdef bool_t allow_fourth_rank

    def __init__(self, int integration_density, integration_volume='octant', bool_t allow_fourth_rank=False):
        """Create the octahedral interpolation averaging scheme for the simulation.

        Args:
            integration_density: The number of triangles along the edge of the octahedron face.
            integration_volume: An enumeration literal, 'octant', 'hemisphere'.
            allow_fourth_rank: Boolean, If True, pre-calculates tables for computing fourth rank tensors.
        """
        self.allow_fourth_rank = allow_fourth_rank
        integration_volume_ = 0
        if integration_volume == 'hemisphere':
            integration_volume_=1
        self.scheme = clib.MRS_create_averaging_scheme(integration_density,
                                    allow_fourth_rank, integration_volume_)

    @property
    def interpolation(self):
        return True

    @property
    def total_orientations(self):
        return self.scheme.total_orientations

    # integration density
    @property
    def integration_density(self):
        return self.scheme.integration_density

    @integration_density.setter
    def integration_density(self, value):
        if isinstance(value, int):
            if value > 0:
                self.scheme = clib.MRS_create_averaging_scheme(value,
                                self.allow_fourth_rank,
                                self.scheme.integration_volume)
                return
        raise ValueError(f"Expecting a positive integer, found {value}.")

    # integration volume
    @property
    def integration_volume(self):
        return __integration_volume_enum_rev__[self.scheme.integration_volume]

    @integration_volume.setter
    def integration_volume(self, value):
        if value in __integration_volume_enum__.keys():
            self.scheme.integration_volume = __integration_volume_enum__[value]
            self.scheme = clib.MRS_create_averaging_scheme(
                    self.scheme.integration_density,
                    self.allow_fourth_rank, value
                )
            return
        raise ValueError(
            (
                "value is not a valid enumeration literal; "
                "permitted: 'octant', 'hemisphere', found {value}.",
            )
        )

    def __reduce_cython__(self):
        return 'AveragingScheme'

    def __eq__(self, other):
        if isinstance(other, AveragingScheme):
            check = [
                self.integration_volume == other.integration_volume,
                self.integration_density == other.integration_density
            ]
            if np.all(check):
                return True
        return False

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return (
            "AveragingScheme(total_orientations={0}, integration_density"
            "={1}, integration_volume={2}, interpolation={3})"
        ).format(
            self.total_orientations,
            self.integration_density,
            self.integration_volume,
            self.interpolation
        )

    def __del__(self):
        clib.MRS_free_averaging_scheme(self.scheme)


cdef class MRSPlan:
    cdef clib.MRS_plan *plan

    def __init__(self, AveragingScheme averaging_scheme, number_of_sidebands,
                rotor_frequency_in_Hz, rotor_angle_in_rad,
                increment, allow_fourth_rank):

        self.plan = clib.MRS_create_plan(averaging_scheme.scheme,
                        number_of_sidebands, rotor_frequency_in_Hz,
                        rotor_angle_in_rad, increment, allow_fourth_rank)

    @property
    def number_of_sidebands(self):
        return self.plan.number_of_sidebands

    @property
    def rotor_frequency_in_Hz(self):
        return self.plan.rotor_frequency_in_Hz

    @property
    def rotor_angle_in_rad(self):
        return self.plan.rotor_angle_in_rad

    # @property
    # def increment(self):
    #     return self.plan.increment

    # @property
    # def allow_fourth_rank(self):
    #     return self.plan.allow_fourth_rank

    # def evaluate(self, R0, R2, R4):
    #     cdef np.ndarray[double complex] R2_c = np.asarray(R2, dtype=np.complex128)
    #     cdef np.ndarray[double complex] R4_c = np.asarray(R4, dtype=np.complex128)
    #     cdef np.ndarray[double complex] output
    #     clib.MRS_get_amplitudes_from_plan(self.plan, &R2_c[0], &R4_c[0])
    #     clib.MRS_get_frequencies_from_plan(self.plan, R0)
    #     # side_band = np.
    #     output = self.plan.vector
    #     return output[::2].reshape(self.plan.number_of_sidebands,
    #                                 self.plan.averaging_scheme.total_orientations)


# @cython.boundscheck(False)
# @cython.wraparound(False)
# def MRS_plan(int integration_density,
#         int number_of_sidebands)
