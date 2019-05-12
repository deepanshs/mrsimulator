
cimport interpolation as clib
cimport numpy as np
import numpy as np
import cython


@cython.boundscheck(False)
@cython.wraparound(False)
def triangle_interpolation(
            freq1,
            freq2,
            freq3,
            offset,
            amp,
            spec,
        ):

        cdef double freq1_c = freq1+0.5
        cdef double freq2_c = freq2+0.5
        cdef double freq3_c = freq3+0.5
        cdef double offset_c = offset
        cdef double amp_c = amp
        cdef np.ndarray[double, ndim=1, mode='c'] spec_c = np.asarray(spec, dtype=np.float64)
        cdef int points_c = spec.size

        clib.triangle_interpolation(&freq1_c,
            &freq2_c,
            &freq3_c,
            &offset_c,
            &amp_c,
            &spec_c[0],
            &points_c,
            )

        return spec_c


