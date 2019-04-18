from libcpp cimport bool, complex

cdef extern from "OCEulerAngle.h":
    ctypedef struct OCEulerAngle:
        double alphaInRadians               # Euler angle alpha
        double betaInRadians                # Euler angle beta
        double gammaInRadians               # Euler angle gamma

    OCEulerAngle OCEulerAngleCreateByAdding(
        OCEulerAngle theAngle1, OCEulerAngle theAngle2
    )

    bool OCEulerAngleEqual(
        OCEulerAngle theAngle1, OCEulerAngle theAngle2
    )

    bool OCEulerAngleZero(OCEulerAngle theAngle)

cdef extern from "MRAngularMomentum.h":
    void full_DLM(double complex *wigner, int l, double * omega)
    double wigner_d(int l, int m1, int m2, double beta)
    double complex DLM(int l, int  m1, int m2, OCEulerAngle omega)

cdef extern from "OCPowderScheme.h":
    ctypedef struct OCPowderScheme:
        OCEulerAngle *angles            # Euler angles
        double *weights                 # Weights of Euler angles
        int size                        # number of Euler angles

    OCPowderScheme OCCreatePowderScheme(int scheme,
                                        int size)