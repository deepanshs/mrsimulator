// -*- coding: utf-8 -*-
//
//  aceraging_scheme.c
//
//  Created by Deepansh J. Srivastava.
//  Contact email = deepansh2012@gmail.com
//

#include "averaging_scheme.h"

/* Free the memory from the mrsimulator plan associated with the spherical
 * averaging scheme */
void MRS_free_averaging_scheme(MRS_averaging_scheme *scheme) {
  free(scheme->amplitudes);
  free(scheme->exp_Im_alpha);
  free(scheme->w2);
  free(scheme->w4);
  free(scheme->wigner_2j_matrices);
  free(scheme->wigner_4j_matrices);
  free(scheme->local_frequency);
  free(scheme->freq_offset);
}

/* Create a new orientation averaging scheme. */
MRS_averaging_scheme *
MRS_create_averaging_scheme(unsigned int geodesic_polyhedron_frequency,
                            bool allow_fourth_rank,
                            unsigned int integration_volume) {

  MRS_averaging_scheme *scheme = malloc(sizeof(MRS_averaging_scheme));

  unsigned int nt = geodesic_polyhedron_frequency;
  unsigned int octant_orientations = (nt + 1) * (nt + 2) / 2;
  unsigned int allocate_size_2, allocate_size_4;

  scheme->octant_orientations = octant_orientations;
  scheme->integration_volume = integration_volume;
  if (integration_volume == 0) {
    scheme->total_orientations = octant_orientations;
  } else if (integration_volume == 1) {
    scheme->total_orientations = 4 * octant_orientations;
  } else {
    scheme->total_orientations = 8 * octant_orientations;
  }

  scheme->geodesic_polyhedron_frequency = nt;

  /* Calculate α, β, and weights over the positive octant. ................. */
  /* ....................................................................... */
  scheme->exp_Im_alpha = malloc_complex128(4 * octant_orientations);
  complex128 *exp_I_beta = malloc_complex128(octant_orientations);
  scheme->amplitudes = malloc_double(octant_orientations);

  octahedron_averaging_setup(geodesic_polyhedron_frequency,
                             &scheme->exp_Im_alpha[3 * octant_orientations],
                             exp_I_beta, scheme->amplitudes);
  /* ----------------------------------------------------------------------- */

  /* calculating exp(-Imα) at every orientation angle α form m=-4 to -1,
   * where α is the azimuthal angle over the positive octant ............... */
  /* ....................................................................... */
  get_exp_Im_alpha(octant_orientations, allow_fourth_rank,
                   scheme->exp_Im_alpha);
  /* ----------------------------------------------------------------------- */

  /**
   * Wigner matrices corresponding to the upper hemisphere.
   *
   * The wigner matrices are evaluated at every β orientation over the
   * positive upper octant.  Note, the β angles from this octant repeat for the
   * other three octant in the upper hemisphere, therfore, only one set of
   * second and fourth rank wigner matrices should suffice.
   */

  // calculating the required space for storing wigner matrices.
  allocate_size_2 = 25 * octant_orientations;
  allocate_size_4 = 81 * octant_orientations;
  if (integration_volume == 2) {
    allocate_size_2 *= 2;
    allocate_size_4 *= 2;
  }

  /* Second rank reduced wigner matrices at every β orientation from the
   * positive upper octant. */
  scheme->wigner_2j_matrices = malloc_double(allocate_size_2);
  wigner_d_matrices_from_exp_I_beta(2, octant_orientations, exp_I_beta,
                                    scheme->wigner_2j_matrices);

  scheme->wigner_4j_matrices = NULL;
  if (allow_fourth_rank) {
    /* Fourth rank reduced wigner matrices at every β orientation from the
     * positive upper octant. */
    scheme->wigner_4j_matrices = malloc_double(allocate_size_4);
    wigner_d_matrices_from_exp_I_beta(4, octant_orientations, exp_I_beta,
                                      scheme->wigner_4j_matrices);
  }

  /**
   * If averaging over a sphere is selected, then calculate the wigner matrices
   * corresponding to the lower hemisphere.
   *
   * The wigner matrices only dependents on the β angles. Going from upper to
   * the lower hemisphere, β -> β+π/2. This implies,
   *      cos(β+π/2) -> -cos(β)
   *      sin(β+π/2) -> sin(β)
   *
   * For evaluating the reduced wigner matrices from the lower hemisphere, the
   * sign of cosine beta is changed. As before, the β angles from any octant
   * from the lower hemisphere repeat for the other three octant in the lower
   * hemisphere, therfore, only one set of second rank and fourth rank reduced
   * wigner matrices should suffice. */
  if (integration_volume == 2) {
    /* cos(beta) is negative in the lower hemisphere */
    cblas_dscal(octant_orientations, -1.0, (double *)exp_I_beta, 2);

    /* Second rank reduced wigner matrices at every β orientation over an octant
     * from the lower hemisphere */
    wigner_d_matrices_from_exp_I_beta(
        2, octant_orientations, exp_I_beta,
        &scheme->wigner_2j_matrices[allocate_size_2]);
    if (allow_fourth_rank) {
      /* Fourth rank reduced wigner matrices at every β orientation. */
      wigner_d_matrices_from_exp_I_beta(
          4, octant_orientations, exp_I_beta,
          &scheme->wigner_4j_matrices[allocate_size_4]);
    }
  }
  free(exp_I_beta);
  /* ----------------------------------------------------------------------- */

  /* Setting up buffers and tables for processing the second rank tensors. . */
  /* ....................................................................... */
  /* w2 is the buffer for storing the frequencies calculated from the
   * second rank tensors. */
  scheme->w2 = malloc_complex128(5 * scheme->total_orientations);

  scheme->w4 = NULL;
  if (allow_fourth_rank) {
    /* w4 is the buffer for storing the frequencies calculated from the
     * fourth rank tensors. */
    scheme->w4 = malloc_complex128(9 * scheme->total_orientations);
  }

  /* buffer to hold the local frequencies and frequency offset. The buffer   *
   * is useful when the rotor angle is off magic angle (54.735 deg). */
  scheme->local_frequency = malloc_double(scheme->total_orientations);
  scheme->freq_offset = malloc_double(octant_orientations);

  return scheme;
}
