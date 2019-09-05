//
//  averaging_scheme.h
//
//  Created by Deepansh J. Srivastava, Sep 3, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#ifndef averaging_scheme_h
#define averaging_scheme_h
#include "config.h"

#include "angular_momentum.h"
#include "powder_setup.h"

/**
 * @struct MRS_averaging_scheme
 * A powder orientation scheme for simulating solid-state NMR line-shapes. An
 * orientation scheme creates buffers and tabulates values for faster
 * computation of bulk NMR line-shapes.
 *
 * The scheme includes
 *   - pre-calculating an array of orientations over the surface of the sphere
 *     where the @f$i^\text{th}@f$ orientation is described by an azimuthal
 *     angle, @f$\alpha_i@f$, a polar angle, @f$\beta_i@f$, and a weighting
 *     factor, @f$w_i@f$,
 *   - pre-calculating arrays of irreducible Wigner small d second-rank,
 *     @f$d^2(n_1, n_2, \beta_i)@f$, and fourth-rank, @f$d^4(n_1, n_2,
 *     \beta_i)@f$, matrices at every orientation angle @f$\beta_i@f$,
 *   - pre-calculating the exponent, @f$\exp(-im\alpha_i)@f$, at every
 *     azimuthal angle, @f$\alpha_i@f$ and for @f$m \in [-4, 0]@f$, and
 *   - allocating buffers for storing and computing frequencies.
 *
 * Creating a new orientation averaging scheme adds an overhead to the
 * computation. Once created, however, the scheme may be re-used for as long as
 * required. This is especially efficient when performing a batch simulation,
 * such as line-shape simulation from thousands of sites.
 */
struct MRS_averaging_scheme {
  unsigned int total_orientations; /**< The total number of orientations. */

  /** \privatesection */
  unsigned short geodesic_polyhedron_frequency; // number of triangles along the
                                                // edge of the octahedron
  unsigned int octant_orientations; //  number of unique orientations on the
                                    //  face of an octant.
  double *amplitudes;         //  array of amplitude scaling per orientation.
  complex128 *exp_Im_alpha;   //  array of cos_alpha per orientation.
  complex128 *w2;             //  buffer for 2nd rank frequency calculation.
  complex128 *w4;             //  buffer for 4nd rank frequency calculation.
  double *wigner_2j_matrices; //  wigner-d 2j matrix per orientation.
  double *wigner_4j_matrices; //  wigner-d 4j matrix per orientation.
  double *local_frequency;    //  buffer for local frequencies.
  double *freq_offset;        //  buffer for local + sideband frequencies.
};

typedef struct MRS_averaging_scheme MRS_averaging_scheme;

/**
 * Create a new orientation averaging scheme.
 *
 * @param geodesic_polyhedron_frequency The value is a positive integer
 *      representing the number of triangles along the edge of an octahedron,
 *      also called the frequency of class I geodesic polyhedra. We use these
 *      polyhedra in calculating the orientation average. Currently, we only
 *      support octahedral as the frequency 1 polyhedra. Higher the geodesic
 *      polyhedron frequency, the closer the polyhedra resemblance a spherical
 *      geometry. For line-shape simulation, a higher geodesic polyhedron
 *      frequency will result in an improved orientation averaging. Read more on
 *      the <a href="https://en.wikipedia.org/wiki/Geodesic_polyhedron">Geodesic
 *      polyhedron</a>.
 *
 * @param allow_fourth_rank If true, the scheme also calculates matrices for
 *            processing fourth rank tensors.
 */
MRS_averaging_scheme *
MRS_create_averaging_scheme(unsigned int geodesic_polyhedron_frequency,
                            bool allow_fourth_rank);

/**
 * Free the memory allocated for the spatial orientation averaging scheme.
 *
 * @param scheme A pointer to the MRS_averaging_scheme.
 */
void MRS_free_averaging_scheme(MRS_averaging_scheme *scheme);

#endif // averaging_scheme_h
