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
 * @param f1 A pointer to the coordinate f11.
 * @param f2 A pointer to the coordinate f12.
 * @param f3 A pointer to the coordinate f13.
 * @param amp A pointer to the area of the vector.
 * @param spec A pointer to the starting index of a one-dimensional array.
 * @param m0 A pointer to the number of points on the 1D grid.
 */
extern void triangle_interpolation1D(double *f1, double *f2, double *f3, double *amp,
                                     double *spec, int *m0);

/**
 * @brief Rasterize a vector triangle with coordinates ((f11, f21), (f12, f22), (f13,
 * f23)) onto a 2D grid.
 *
 * @param f11 A pointer to the coordinate f11.
 * @param f12 A pointer to the coordinate f12.
 * @param f13 A pointer to the coordinate f13.
 * @param f21 A pointer to the coordinate f21.
 * @param f22 A pointer to the coordinate f22.
 * @param f23 A pointer to the coordinate f23.
 * @param amp A pointer to the area of the vector.
 * @param spec A pointer to the starting index of a two-dimensional array.
 * @param m0 An interger with the rows in the 2D grid.
 * @param m1 An interger with the columns in the 2D grid.
 */
extern int triangle_interpolation2D(double *f11, double *f12, double *f13, double *f21,
                                    double *f22, double *f23, double *amp, double *spec,
                                    int m0, int m1);

/**
 * @brief Sum amplitudes from the triangles interpolations over the region of an octant.
 * The samplings over the octant is as per Alderman and Grand scheme.
 *
 * @param nt Number of triangles along the edge of the octant.
 * @param freq A pointer to an array of frequencies evaluated at octant coordinates.
 * @param amp A pointer to the amplitudes for the frequencies at octant coordinates.
 * @param stride Stride setp for the amplitudes (amp) array.
 * @param n_spec Number of points in the spectrum array (spec)
 * @param spec A pointer to the starting index of a one-dimensional array
 */
void octahedronDeltaInterpolation(const unsigned int nt, double *freq, double *amp,
                                  int stride, int n_spec, double *spec);

extern void octahedronInterpolation(double *spec, double *freq, const unsigned int nt,
                                    double *amp, int stride, int m);

extern void octahedronInterpolation2D(double *spec, double *freq1, double *freq2,
                                      int nt, double *amp, int stride, int m0, int m1);
