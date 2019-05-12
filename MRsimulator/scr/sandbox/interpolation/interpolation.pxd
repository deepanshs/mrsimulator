cdef extern from "powder_setup.h":
    int triangle_interpolation(double *freq1,
            double *freq2,
            double *freq3,
            double *offset,
            double *amp,
            double *spec,
            int *points
            )
            