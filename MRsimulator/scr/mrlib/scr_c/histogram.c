#include "histogram.h"
#include <math.h>

void histogram1d(
					double * hist_amp,

					double * x,
					double * weights,
					int size_x,

					// freq sprectum information
					double freq_min,
					double freq_max,
					int size_hist // also equal to freq_size
				){

	// The function evaluates the histogram of values in x along
	// with the rescpective weights in `weights`, on a scale that
	// is bounded by `xmin` and `xmax`, and discretized into
	// m points. the scale is hist.

	int i, ix;
	double tx, normx;
	double sampling_interval = (freq_max - freq_min) / (double)(size_hist - 1.0);
	double half_sampling_interval = 0.5 * sampling_interval;

	normx = 1.0/sampling_interval;

	for(i = 0; i<=size_x-1; i++){
		tx = x[i];
		if (tx >= freq_min - half_sampling_interval && \
				tx < freq_max + half_sampling_interval ) {
					ix = floor((tx - freq_min) * normx + 0.5);
					hist_amp[ix] += weights[i];
		}
	}
}





// void histogram2d(double * hist,
//                  double * x,
//                  double * y,
//                  double * weights,
//                  double xmin,
//                  double xmax,
//                  int nx,
//                  double ymin,
//                  double ymax,
//                  int ny
//                  int m)

// subroutine  histogram2d(hist, x, y, weight, xmin, xmax, nx, &
//                         ymax, ymin, ny, m, npros)

// use omp_lib
// implicit none

// integer*4, intent(in)    :: nx, ny, m, npros
// !f2py integer*4, intent(in) :: nx, ny, m, npros

// double precision, dimension(0:m-1), intent(in)   :: x, y, weight
// !f2py double precision, intent(in) :: x, y, weight

// double precision, intent(in)   :: xmin, xmax, ymin, ymax
// !f2py double precision, intent(in) :: xmin, xmax, ymin, ymax


// double precision, dimension(0:nx*ny-1), intent(out) :: hist
// !f2py double precision, intent(out) :: hist

// integer*4 :: i, ix, iy, index
// double precision :: tx, ty, normx, normy, sampling_interval_x, sampling_interval_y

// double *hist = create1DDoubleArray(nx*ny-1)

// sampling_interval_x = (xmax - xmin) / (nx - 1)
// sampling_interval_y = (ymax - ymin) / (ny - 1)

// normx = 1.0/sampling_interval_x
// normy = 1.0/sampling_interval_y

// !x = int(x*normx - xmin*normx)
// !call

// for(i = 0; i<=m-1; i++){
//     tx = x[i]
//     ty = y[i]
//     if (tx >= xmin - 0.5 * sampling_interval_x && \
//         tx < xmax + 0.5 * sampling_interval_x  && \
//         ty >= ymin - 0.5 * sampling_interval_y && \
//         ty < ymax + 0.5 * sampling_interval_y ) {
//             ix = (int)((tx - xmin) * normx + 0.5 * sampling_interval_x)
//             iy = (int)((ty - ymin) * normy + 0.5 * sampling_interval_y)
//             index = iy + ny * ix
//             hist(index) = hist(index) + weight(i)
//     }
// }
// !$omp end do
// !$omp end parallel

// end subroutine  histogram2d
