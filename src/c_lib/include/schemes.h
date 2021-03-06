// -*- coding: utf-8 -*-
//
//  schemes.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Sep 3, 2019.
//  Contact email = srivastava.89@osu.edu
//

#ifndef averaging_scheme_h
#define averaging_scheme_h
#include "angular_momentum.h"
#include "config.h"
#include "octahedron.h"

/**
 * @struct MRS_averaging_scheme
 * A powder orientation scheme for simulating solid-state NMR spectra. An orientation
 * scheme creates buffers and tabulates values for faster computation of bulk NMR
 * spectra.
 *
 * The scheme includes
 *   - pre-calculating an array of orientations over the surface of the sphere where the
 *     @f$i^\text{th}@f$ orientation is described by an azimuthal angle, @f$\alpha_i@f$,
 * a polar angle, @f$\beta_i@f$, and a weighting factor, @f$w_i@f$,
 *   - pre-calculating arrays of irreducible Wigner small d second-rank, @f$d^2(n_1,
 * n_2, \beta_i)@f$, and fourth-rank, @f$d^4(n_1, n_2, \beta_i)@f$, matrices at every
 * orientation angle @f$\beta_i@f$,
 *   - pre-calculating the exponent, @f$\exp(-im\alpha_i)@f$, at every azimuthal angle,
 *     @f$\alpha_i@f$ and for @f$m \in [-4, 0]@f$, and
 *   - allocating buffers for storing and computing frequencies.
 *
 * Creating a new orientation averaging scheme adds an overhead to the computation. Once
 * created, however, the scheme may be re-used for as long as required. This is
 * especially efficient when performing a batch simulation, such as simulations from
 * thousands of sites.
 */
typedef struct MRS_averaging_scheme {
  unsigned int total_orientations; /**< The total number of orientations. */

  /** \privatesection */
  unsigned int integration_density;  //  # triangles along the edge of the octahedron.
  unsigned int integration_volume;   //  0-octant, 1-hemisphere, 2-sphere.
  unsigned int octant_orientations;  //  # unique orientations on the face of an octant.
  double *amplitudes;                //  array of amplitude scaling per orientation.
  complex128 *exp_Im_alpha;          //  array of cos_alpha per orientation.
  complex128 *w2;                    //  buffer for 2nd rank frequency calculation.
  complex128 *w4;                    //  buffer for 4nd rank frequency calculation.
  double *wigner_2j_matrices;        //  wigner-d 2j matrix per orientation.
  double *wigner_4j_matrices;        //  wigner-d 4j matrix per orientation.
  bool allow_fourth_rank;  //  If true, compute wigner matrices for wigner-d 4j.
} MRS_averaging_scheme;

// typedef struct MRS_averaging_scheme;

/**
 * Create a new orientation averaging scheme.
 *
 * @param integration_density The value is a positive integer representing the number of
 * triangles along the edge of an octahedron, also called the frequency of class I
 * geodesic polyhedra. We use these polyhedra in calculating the orientation average.
 * Currently, we only support octahedral as the frequency 1 polyhedra. Higher the
 * geodesic polyhedron frequency, the closer the polyhedra resemblance a spherical
 * geometry. For spectrum simulation, a higher geodesic polyhedron frequency will result
 * in an improved orientation averaging. Read more on the <a
 * href="https://en.wikipedia.org/wiki/Geodesic_polyhedron">Geodesic polyhedron</a>.
 *
 * @param allow_fourth_rank If true, the scheme also calculates matrices for fourth-rank
 * tensors.
 * @param integration_volume An enumeration. 0=octant, 1=hemisphere
 */
MRS_averaging_scheme *MRS_create_averaging_scheme(unsigned int integration_density,
                                                  bool allow_fourth_rank,
                                                  unsigned int integration_volume);

/**
 * Create a new orientation averaging scheme from given alpha and beta.
 *
 * @param alpha A pointer to an array of size `n_angles` holding the alpha values of
 * type double.
 * @param beta A pointer to an array of size `n_angles` holding the beta values of type
 * double.
 * @param weight A pointer to an array of size `n_angles` holding the weights of the
 * angle pair (alpha, beta) of type double.
 * @param n_angles An unsigned int equal to the total number of alpha, beta, and weight
 * values.
 * @param allow_fourth_rank If true, the scheme also calculates matrices for fourth-rank
 * tensors.
 */
MRS_averaging_scheme *MRS_create_averaging_scheme_from_alpha_beta(
    double *alpha, double *beta, double *weight, unsigned int n_angles,
    bool allow_fourth_rank);

/**
 * Free the memory allocated for the spatial orientation averaging scheme.
 *
 * @param scheme A pointer to the MRS_averaging_scheme.
 */
void MRS_free_averaging_scheme(MRS_averaging_scheme *scheme);

#endif  // averaging_scheme_h

#ifndef fftw_scheme_h
#define fftw_scheme_h
#include "fftw3.h"

typedef struct MRS_fftw_scheme {
  /** \privatesection */
  /** The buffer to hold the sideband amplitudes as stride 2 array after mrsimulator
   * processing. */
  fftw_complex *vector;     // holds the amplitude of sidebands.
  fftw_plan the_fftw_plan;  //  The plan for fftw routine.
} MRS_fftw_scheme;

MRS_fftw_scheme *create_fftw_scheme(unsigned int total_orientations,
                                    unsigned int number_of_sidebands);

void MRS_free_fftw_scheme(MRS_fftw_scheme *fftw_scheme);

#endif  // fftw_scheme_h
