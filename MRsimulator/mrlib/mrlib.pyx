from mrlib cimport (
    # Euler angle
    OCEulerAngle,
    OCEulerAngleCreateByAdding,

    # angular momentum
    full_DLM,
    wigner_d,
    DLM,

    # powder scheme
    OCPowderScheme
)

from unit import string_to_quantity

from cpython cimport array as array
cimport numpy as np
import numpy as np
from libc.stdlib cimport free

__doc__  =  """
An NMR lineshape simulator module.
The module supports the following:

    - single spin static chemical shift anisotropy (CSA) lineshape.
    - single spin chemical shift anisotropy (CSA) spinning sideband.
"""


cdef class EulerAngles:
    r"""
    Create an instance of an EulerAngles class.

    The three Euler angles are alpha, beta and gamma.

    :params: alpha: string: alpha with unit.
    :params: beta: string: beta with unit.
    :params: gamma: string: gamma with unit.

    :returns: An instance of ``EulerAngles``.
    """

    cdef OCEulerAngle c_angle

    cdef dict _angle

    def __cinit__(
            self,
            alpha='0',
            beta='0',
            gamma='0'
            ):
        
        self._angle = {
            'alpha': string_to_quantity(alpha).to('rad'),
            'beta': string_to_quantity(beta).to('rad'),
            'gamma': string_to_quantity(gamma).to('rad')
        }

        alpha = self._angle['alpha'].value
        beta = self._angle['beta'].value
        gamma = self._angle['gamma'].value

        self.c_angle = OCEulerAngle(
                            alpha,
                            beta,
                            gamma
                        )

    @staticmethod
    cdef create(OCEulerAngle angle):
        new_angle = EulerAngles(
                        alpha=str(angle.alphaInRadians) + ' rad',
                        beta=str(angle.betaInRadians) + ' rad',
                        gamma=str(angle.gammaInRadians) + ' rad'
                    )
        return new_angle

    # Define the properties of the c structure
# alpha
    @property
    def alpha(self):
        """Return the alpha value of the EulerAngles instance."""
        return self._angle['alpha']

    @alpha.setter
    def alpha(self, value):
        self.c_angle.alphaInRadians = string_to_quantity(value).to('rad').value

# beta
    @property
    def beta(self):
        """Return the beta value of the EulerAngles instance."""
        return self._angle['beta']

    @beta.setter
    def beta(self, value):
        self.c_angle.betaInRadians = string_to_quantity(value).to('rad').value

# gamma
    @property
    def gamma(self):
        return self._angle['gamma']

    @gamma.setter
    def gamma(self, value):
        """Return the gamma value of the EulerAngles instance."""
        self.c_angle.gammaInRadians = string_to_quantity(value).to('rad').value 

# methods 
# print 
    def __str__(self):
        return str(self._angle)

# add two Euler angles
    def __add__(EulerAngles self, EulerAngles other):
        new_c_angle = OCEulerAngleCreateByAdding(
                            self.c_angle,
                            other.c_angle
                        )
        return EulerAngles.create(new_c_angle)


cdef class EulerAnglesInRadians:
    """
    Create an instance of an EulerAnglesInRadians.

    The three Euler angles are alpha, beta and gamma.

    :params: alpha: float: alpha in radians.
    :params: beta: float: beta in radians.
    :params: gamma: float: gamma in radians.

    :returns: An instance of ``EulerAnglesInRadians``.
    """

    cdef OCEulerAngle c_angle

    def __cinit__(
            self,
            alpha=0,
            beta=0,
            gamma=0
            ):

        self.c_angle = OCEulerAngle(
                            alpha,
                            beta,
                            gamma
                        )

    @staticmethod
    cdef create(OCEulerAngle angle):
        new_angle = EulerAnglesInRadians()
        new_angle.c_angle = angle
        return new_angle

    # Define the properties of the c structure
# alpha
    @property
    def alpha(self):
        """Return the alpha value of the EulerAnglesInRadians instance."""
        return self.c_angle.alphaInRadians

    @alpha.setter
    def alpha(self, value):
        self.c_angle.alphaInRadians = value

# beta
    @property
    def beta(self):
        """Return the beta value of the EulerAnglesInRadians instance."""
        return self.c_angle.betaInRadians

    @beta.setter
    def beta(self, value):
        self.c_angle.betaInRadians = value

# gamma
    @property
    def gamma(self):
        return self.c_angle.alphaInRadians

    @gamma.setter
    def gamma(self, value):
        """Return the gamma value of the EulerAnglesInRadians instance."""
        self.c_angle.gammaInRadians = value

# methods 
# print 
    def __str__(self):
        return str(self.c_angle)

# add two Euler angles
    def __add__(EulerAnglesInRadians self, EulerAnglesInRadians other):
        new_c_angle = OCEulerAngleCreateByAdding(
                            self.c_angle,
                            other.c_angle
                        )
        return EulerAnglesInRadians.create(new_c_angle)


# cdef class WrappingList:
#     cdef int size
#     cdef OCEulerAngle *angles
#     cdef double *weights

#     def __cinit__(self):
#         self.size=0
#         self.angles=NULL
#         self.weights=NULL

#     def __dealloc__(self):
#         free(self.angles)
#         free(self.weights)
#         print "deallocated"#just a check

#     def __getitem__(self, index):
#         if index<0 or index>=self.size:
#             raise IndexError("list index out of range")
#         return EulerAnglesInRadians(self.angles[index].alphaInRadians,
#                                     self.angles[index].betaInRadians,
#                                     self.angles[index].gammaInRadians)


# # def create_structs(token, size):
# #     cdef OCPowderScheme res = OCCreatePowderScheme(token, size)
# #     lst=WrappingList()
# #     lst.size, lst.arr=res.size, res.arr
# #     return lst 

# cdef class PowderScheme:
#     """
#     Create an instance of an PowderScheme.

#     The are following powder schemes
#         1) lebedev

#     :attr: scheme: string: A valid name of the scheme.
#     :attr: size: int: The number of Euler Angles generated from the scheme.
#     :attr: angles: EulerAnglesInRadians: The Euler angles.

#     :returns: An instance of ``PowderScheme``.
#     """

#     cdef OCPowderScheme *ptr   # pointer to powder scheme
#     # cdef double *ptr_weights        # pointer to weights

#     def __init__(self):
#        self.ptr = NULL

#     @property
#     def angles(self):
#         return EulerAnglesInRadians(self.ptr.angles.alphaInRadians,
#                                     self.ptr.angles.betaInRadians,
#                                     self.ptr.angles.gammaInRadians)

#     # @mult.setter
#     # def mult(self, value):
#     #     self.ptr.mult = value

#     @property
#     def weigths(self):
#         return self.ptr.weights

#     @add.setter
#     def add(self, value):
#         self.ptr.add = value
#     # cdef OCPowderScheme res = OCCreatePowderScheme(token, size)
#     # lst=WrappingList()
#     # lst.size, lst.arr=res.size, res.arr
#     # return lst 

    
#     # cdef OCEulerAngle *c_angle
#     # cdef int size

#     def __cinit__(
#             self,
#             scheme='',
#             size=48,
#             ):

#         if scheme == 'lebedev': token = 1

#         cdef OCPowderScheme c_scheme = OCCreatePowderScheme(token, size)
#         self.ptr = *c_scheme
#         self.lst = WrappingList()
#         self.lst.size = c_scheme.size
#         self.lst.angles = c_scheme.angles
#         self.lst.weights = c_scheme.weights
#         # self.size = size

#     @property
#     def angles(self):
#         return self.lst.angles
#         # EulerAnglesInRadians angle[self.size]
#         # cdef array.array amp = array.array(EulerAnglesInRadians, np.empty(self.size))
#         # # EulerAnglesInRadians angle[self.size]
#         # angle.c_angle = self.c_scheme.angles
#         # return 

#     @property
#     def weights(self):
#         return self.lst.weights
#         # cdef array.array amp = array.array('d', np.empty(self.size))
#         # amp.data.as_doubles = self.c_scheme.weights
#         # print (amp)
#         # return np.asarray(amp)




def wigner_D_matrix(int l=2, EulerAnglesInRadians omega=EulerAnglesInRadians(0.,0.,0.)):
    r"""Return the Wigner :math:`d^l(\beta)` matrix."""
    cdef int size = 2*l+1
    cdef np.ndarray[complex, ndim=1, mode='c'] matrix = np.zeros(size**2).astype(np.complex128)
    cdef np.ndarray[double, ndim=1, mode='c'] om = np.asarray([omega.alphaInRadians, omega.betaInRadians, omega.gammaInRadians])

    full_DLM(
        &matrix[0],
        l,
        &om[0]
    )
    return matrix.reshape(size, size)


def wigner_d_element(int l, int m1, int m2, double beta=0):
    r"""Return the Wigner :math:`d^l_{m_1,m_2}(\beta)` element."""

    return wigner_d(l, m1, m2, beta)


def wigner_D_element(int l, int m1, int m2, EulerAnglesInRadians omega):
    r"""Return the Wigner :math:`D^l_{m_1,m_2}(\alpha, \beta, \gamma)` element."""

    return DLM(l, m1, m2, omega.c_angle)
