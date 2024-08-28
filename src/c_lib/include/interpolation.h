// -*- coding: utf-8 -*-
//
//  interpolation.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Apr 11, 2019.
//  Contact email = srivastava.89@osu.edu
//

#include "config.h"

/**
 * @brief Create a triangle with coordinates (f1, f2, f2) onto a 1D grid.
 *
 * @param f1 Coordinate f11.
 * @param f2 Coordinate f12.
 * @param f3 Coordinate f13.
 * @param amp Area of triangle.
 * @param spec A pointer to the starting index of a one-dimensional array.
 * @param m0 The number of points on the 1D grid.
 * @param type The type of interpolation for the isotropic components.
 *          0. delta-interpolation,
 *          1. gaussian-interpolation.
 */
extern void triangle_interpolation1D(double f1, double f2, double f3, double amp,
                                     double *spec, int m0, unsigned int iso_intrp);

extern void triangle_interpolation1D_linear(double f1, double f2, double f3, double amp,
                                            double *spec, int m0);

extern void one_d_averaging(double *spec, const unsigned int freq_size, double *freq,
                            double *amp_real, double *amp_imag, int dimension_count,
                            const unsigned int position_size, int32_t *positions,
                            const unsigned int nt, bool user_defined,
                            bool interpolation);

extern void two_d_averaging(double *spec, const unsigned int freq_size, double *freq1,
                            double *freq2, double *amp, int amp_stride,
                            const unsigned int position_size, int32_t *positions,
                            int dimension0_count, int dimension1_count,
                            unsigned int iso_intrp, const unsigned int nt,
                            bool user_defined, bool interpolation);

extern void hist1d(double *spec, const unsigned int freq_size, double *freq,
                   double *amp, int m);

extern void hist2d(double *spec, const unsigned int freq_size, double *freq_1,
                   double *freq_2, double *amp, int amp_stride, int m0, int m1);

extern void generic_2d_triangle_average(double *spec, const unsigned int freq_size,
                                        double *freq1, double *freq2, double *amp,
                                        int amp_stride, int m0, int m1,
                                        const unsigned int position_size,
                                        int32_t *positions, unsigned int iso_intrp);

extern void generic_1d_triangle_average(double *spec, const unsigned int freq_size,
                                        double *freq, double *amp, int m,
                                        const unsigned int position_size,
                                        int32_t *positions);

extern void triangle_interpolation1D_gaussian(double f1, double f2, double f3,
                                              double amp, double *spec, int m0);
/**
 * @brief Rasterize a vector triangle with coordinates ((f11, f21), (f12, f22), (f13,
 * f23)) onto a 2D grid.
 *
 * @param f11 Coordinate f11.
 * @param f12 Coordinate f12.
 * @param f13 Coordinate f13.
 * @param f21 Coordinate f21.
 * @param f22 Coordinate f22.
 * @param f23 Coordinate f23.
 * @param amp Area of 2D triangle.
 * @param spec A pointer to the starting index of a two-dimensional array.
 * @param m0 An integer with the rows in the 2D grid.
 * @param m1 An integer with the columns in the 2D grid.
 * @param iso_intrp Linear=0 | Gaussian=1 isotropic interpolation scheme.
 */
extern void triangle_interpolation2D(double f11, double f12, double f13, double f21,
                                     double f22, double f23, double amp, double *spec,
                                     int m0, int m1, unsigned int iso_intrp);

/**
 * @brief Sum amplitudes from the triangles interpolations over the region of an octant.
 * The samplings over the octant is as per Alderman and Grand scheme.
 *
 * @param nt Number of triangles along the edge of the octant.
 * @param freq Delta frequency.
 * @param amp A pointer to the amplitudes for the frequencies at octant coordinates.
 * @param stride Stride setp for the amplitudes (amp) array.
 * @param n_spec Number of points in the spectrum array (spec)
 * @param spec A pointer to the starting index of a one-dimensional array
 * @param iso_intrp Linear=0 | Gaussian=1 isotropic interpolation scheme.
 */
void octahedronDeltaInterpolation(const unsigned int nt, double freq, double *amp,
                                  const unsigned int stride, int n_spec, double *spec,
                                  unsigned int iso_intrp);

extern void octahedronInterpolation(double *spec, double *freq, const unsigned int nt,
                                    double *amp, const unsigned int stride, int m);

extern void octahedronInterpolation2D(double *spec, double *freq1, double *freq2,
                                      const unsigned int nt, double *amp,
                                      const unsigned int stride, int m0, int m1,
                                      unsigned int iso_intrp);
