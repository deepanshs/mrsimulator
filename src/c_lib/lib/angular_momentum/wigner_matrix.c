// -*- coding: utf-8 -*-
//
//  wigner_matrix.c
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Mar 9, 2021.
//  Contact email = srivastava.89@osu.edu
//

#include "angular_momentum/wigner_matrix.h"

complex128 IOTA = {0.0, 1.0};
complex128 NEGATIVE_IOTA = {0.0, -1.0};

/* calculate Wigner rotation matrices */

// ✅ .. note: (wigner_d_matrices) tested with pytest
// .........................
void wigner_d_matrices(const int l, const int n, const double *beta, double *wigner) {
  complex128 *exp_I_beta = malloc_complex128(n);
  vm_cosine_I_sine(n, beta, exp_I_beta);
  wigner_d_matrices_from_exp_I_beta(l, n, false, exp_I_beta, wigner);
  free(exp_I_beta);
}

// Compute full wigner 2j (5 x 5) matrix.
static inline void wigner_2j(const int n, const void *restrict exp_I_beta,
                             double *wigner) {
  double cx, cx2, sx, temp, *exp_I_beta_ = (double *)exp_I_beta, t1;
  int i;
  for (i = 0; i < n; i++) {
    cx = *exp_I_beta_++;
    sx = *exp_I_beta_++;
    cx2 = cx * cx;

    t1 = (1. + cx);
    temp = -sx * t1 * 0.5;
    wigner[19] = temp;   //  2,  1 // 19
    wigner[5] = -temp;   // -2, -1 //  5
    wigner[23] = -temp;  //  1,  2 // 23
    wigner[1] = temp;    // -1, -2 //  1

    temp = t1 * t1 * 0.25;
    wigner[24] = temp;  //  2,  2 // 24
    wigner[0] = temp;   // -2, -2 //  0

    t1 = (1. - cx);
    temp = -sx * t1 * 0.5;
    wigner[9] = temp;    //  2, -1 //  9
    wigner[15] = -temp;  // -2,  1 // 15
    wigner[3] = temp;    //  1, -2 //  3
    wigner[21] = -temp;  // -1,  2 // 21

    temp = t1 * t1 * 0.25;
    wigner[4] = temp;   //  2, -2 //  4
    wigner[20] = temp;  // -2,  2 // 20

    temp = 0.6123724355 * sx * sx;
    wigner[14] = temp;  //  2,  0 // 14
    wigner[10] = temp;  // -2,  0 // 10
    wigner[22] = temp;  //  0,  2 // 22
    wigner[2] = temp;   //  0, -2 //  2

    temp = 1.224744871 * sx * cx;
    wigner[13] = -temp;  //  1,  0 // 13
    wigner[17] = temp;   //  0,  1 // 17
    wigner[7] = -temp;   //  0, -1 //  7
    wigner[11] = temp;   // -1,  0 // 11

    temp = (2.0 * cx2 + cx - 1.) * 0.5;
    wigner[18] = temp;  //  1,  1 // 18
    wigner[6] = temp;   // -1, -1 //  6

    temp = -(2.0 * cx2 - cx - 1.) * 0.5;
    wigner[8] = temp;   //  1, -1 //  8
    wigner[16] = temp;  // -1,  1 // 16

    wigner[12] = 1.5 * cx2 - 0.5;  // 0,  0 // 12

    wigner += 25;  // increment counter to next matrix
  }
}

// Compute half wigner 4j (5 x 3) matrix.
static inline void wigner_2j_half(const int n, const void *restrict exp_I_beta,
                                  double *wigner) {
  double cx, cx2, sx, temp, *exp_I_beta_ = (double *)exp_I_beta, t1;
  int i;
  for (i = 0; i < n; i++) {
    cx = *exp_I_beta_++;
    sx = *exp_I_beta_++;
    cx2 = cx * cx;

    t1 = (1. + cx);
    temp = t1 * t1 * 0.25;
    wigner[0] = temp;  // -2, -2 //  0

    temp = -sx * t1 * 0.5;
    wigner[1] = temp;   // -1, -2 //  1
    wigner[5] = -temp;  // -2, -1 //  5

    t1 = (1. - cx);
    temp = -sx * t1 * 0.5;
    wigner[3] = temp;  //  1, -2 //  3
    wigner[9] = temp;  //  2, -1 //  9

    temp = t1 * t1 * 0.25;
    wigner[4] = temp;  //  2, -2 //  4

    temp = 0.6123724355 * sx * sx;
    wigner[2] = temp;   //  0, -2 //  2
    wigner[14] = temp;  //  2,  0 // 14
    wigner[10] = temp;  // -2,  0 // 10

    temp = 1.224744871 * sx * cx;
    wigner[7] = -temp;   //  0, -1 //  7
    wigner[11] = temp;   // -1,  0 // 11
    wigner[13] = -temp;  //  1,  0 // 13

    temp = (2.0 * cx2 + cx - 1.) * 0.5;
    wigner[6] = temp;  // -1, -1 //  6

    temp = -(2.0 * cx2 - cx - 1.) * 0.5;
    wigner[8] = temp;  //  1, -1 //  8

    wigner[12] = 1.5 * cx2 - 0.5;  // 0,  0 // 12

    wigner += 15;  // increment counter to next matrix
  }
}

// Compute full wigner 4j (9 x 9) matrix.
static inline void wigner_4j(const int n, const void *restrict exp_I_beta,
                             double *wigner) {
  double cx, cx2, sx, temp, *exp_I_beta_ = (double *)exp_I_beta;
  double sx2, sx3, cxp1, cxm1, cxp12, cxm12, cxm13, cxp13;
  int i;
  for (i = 0; i < n; i++) {
    cx = *exp_I_beta_++;
    sx = *exp_I_beta_++;
    cx2 = cx * cx;

    sx2 = sx * sx;
    sx3 = sx2 * sx;

    cxp1 = (1. + cx);
    cxm1 = (1. - cx);
    cxp12 = cxp1 * cxp1;
    cxm12 = cxm1 * cxm1;
    cxm13 = cxm12 * cxm1;
    cxp13 = cxp12 * cxp1;

    temp = 0.0625 * cxp12 * cxp12;
    wigner[0] = temp;   // -4, -4 //  0
    wigner[80] = temp;  //  4,  4 // 80

    temp = 0.0625 * cxm12 * cxm12;
    wigner[72] = temp;  // -4,  4 // 72
    wigner[8] = temp;   //  4, -4 //  8

    temp = -0.1767766953 * cxp13 * sx;
    wigner[1] = temp;    // -3, -4 //  1
    wigner[9] = -temp;   // -4, -3 //  9
    wigner[79] = -temp;  //  3,  4 // 79
    wigner[71] = temp;   //  4,  3 //  9

    temp = -0.1767766953 * cxm13 * sx;
    wigner[7] = temp;    //  3, -4 //  7
    wigner[63] = -temp;  // -4,  3 // 63
    wigner[73] = -temp;  // -3,  4 // 73
    wigner[17] = temp;   //  4, -3 // 17

    temp = -0.4677071733 * cxp1 * sx3;
    wigner[53] = temp;   //  4,  1 // 53
    wigner[27] = -temp;  // -4, -1 // 27
    wigner[77] = -temp;  //  1,  4 // 77
    wigner[3] = temp;    // -1, -4 //  3

    temp = -0.4677071733 * cxm1 * sx3;
    wigner[35] = temp;   //  4, -1 // 35
    wigner[45] = -temp;  // -4,  1 // 45
    wigner[75] = -temp;  // -1,  4 // 75
    wigner[5] = temp;    //  1, -4 //  5

    temp = 0.5229125166 * sx3 * sx;
    wigner[44] = temp;  //  4,  0 // 44
    wigner[36] = temp;  // -4,  0 // 36
    wigner[76] = temp;  //  0,  4 // 76
    wigner[4] = temp;   //  0, -4 //  4

    temp = -1.4790199458 * sx3 * cx;
    wigner[43] = temp;   //  3,  0 // 43
    wigner[37] = -temp;  // -3,  0 // 37
    wigner[67] = -temp;  //  0,  3 // 67
    wigner[13] = temp;   //  0, -3 // 13

    temp = 0.3307189139 * sx2 * cxp12;
    wigner[78] = temp;  //  2,  4 // 78
    wigner[2] = temp;   // -2, -4 //  2
    wigner[62] = temp;  //  4,  2 // 62
    wigner[18] = temp;  // -4, -2 // 18

    temp = 0.3307189139 * sx2 * cxm12;
    wigner[6] = temp;   //  2, -4 //  6
    wigner[74] = temp;  // -2,  4 // 74
    wigner[54] = temp;  // -4,  2 // 54
    wigner[26] = temp;  //  4, -2 // 26

    temp = 0.4677071733 * cxp12 * sx * (2. * cx - 1.);
    wigner[69] = temp;   //  2,  3 // 69
    wigner[11] = -temp;  // -2, -3 // 11
    wigner[61] = -temp;  //  3,  2 // 61
    wigner[19] = temp;   // -3, -2 // 19

    temp = 0.4677071733 * cxm12 * sx * (-2. * cx - 1.);
    wigner[15] = temp;   //  2, -3 // 15
    wigner[65] = -temp;  // -2,  3 // 65
    wigner[55] = -temp;  // -3,  2 // 55
    wigner[25] = temp;   //  3, -2 // 25

    temp = 0.25 * cxp12 * (1. - 7. * cxm1 + 7. * cxm12);
    wigner[60] = temp;  //  2,  2 // 60
    wigner[20] = temp;  // -2, -2 // 20

    temp = 0.25 * cxm12 * (1. - 7. * cxp1 + 7. * cxp12);
    wigner[56] = temp;  // -2,  2 // 56
    wigner[24] = temp;  //  2, -2 // 24

    temp = 0.3952847075 * sx2 * (7. * cx2 - 1.0);
    wigner[42] = temp;  //  2,  0 // 42
    wigner[38] = temp;  // -2,  0 // 38
    wigner[58] = temp;  //  0,  2 // 58
    wigner[22] = temp;  //  0, -2 // 22

    temp = 0.125 * cxp13 * (-3. + 4. * cx);
    wigner[10] = temp;  // -3, -3 // 10
    wigner[70] = temp;  //  3,  3 // 70

    temp = 0.125 * cxm13 * (3. + 4. * cx);
    wigner[64] = temp;  // -3,  3 // 64
    wigner[16] = temp;  //  3, -3 // 16

    temp = 0.3307189139 * cxm1 * cxp12 * (-1. + 4. * cx);
    wigner[12] = temp;  // -1, -3 // 12
    wigner[28] = temp;  // -3, -1 // 28
    wigner[68] = temp;  //  1,  3 // 68
    wigner[52] = temp;  //  3,  1 // 52

    temp = 0.3307189139 * cxm12 * cxp1 * (1. + 4. * cx);
    wigner[14] = temp;  //  1, -3 // 14
    wigner[46] = temp;  // -3,  1 // 46
    wigner[66] = temp;  // -1,  3 // 66
    wigner[34] = temp;  //  3, -1 // 34

    temp = -0.5590169944 * (4. - 18. * cxm1 + 21. * cxm12 - 7. * cxm13) * sx;
    wigner[41] = temp;   //  1,  0 // 41
    wigner[39] = -temp;  // -1,  0 // 39
    wigner[49] = -temp;  //  0,  1 // 49
    wigner[31] = temp;   //  0, -1 // 31

    temp = -0.3535533906 * (3. - 10.5 * cxm1 + 7. * cxm12) * sx * cxp1;
    wigner[51] = temp;   //  2,  1 // 51
    wigner[29] = -temp;  // -2, -1 // 29
    wigner[59] = -temp;  //  1,  2 // 59
    wigner[21] = temp;   // -1, -2 // 21

    temp = -0.3535533906 * (10. - 17.5 * cxm1 + 7. * cxm12) * sx * cxm1;
    wigner[23] = temp;   //  1, -2 // 23
    wigner[57] = -temp;  // -1,  2 // 57
    wigner[47] = -temp;  // -2,  1 // 47
    wigner[33] = temp;   //  2, -1 // 33

    temp = 0.5 * (1. - 9. * cxm1 + 15.75 * cxm12 - 7. * cxm13) * cxp1;
    wigner[30] = temp;  // -1, -1 // 30
    wigner[50] = temp;  //  1,  1 // 50

    temp = 0.5 * (10. - 30. * cxm1 + 26.25 * cxm12 - 7. * cxm13) * cxm1;
    wigner[32] = temp;  //  1, -1 // 32
    wigner[48] = temp;  // -1,  1 // 48

    temp = 0.125 * (3. - 30. * cx2 + 35.0 * cx2 * cx2);
    wigner[40] = temp;  //  0,  0 // 40

    wigner += 81;  // increment counter to next matrix
  }
}

// Compute half wigner 4j (9 x 5) matrix.
static inline void wigner_4j_half(const int n, const void *restrict exp_I_beta,
                                  double *wigner) {
  double cx, cx2, sx, temp, *exp_I_beta_ = (double *)exp_I_beta;
  double sx2, sx3, cxp1, cxm1, cxp12, cxm12, cxm13, cxp13;
  int i;
  for (i = 0; i < n; i++) {
    cx = *exp_I_beta_++;
    sx = *exp_I_beta_++;
    cx2 = cx * cx;

    sx2 = sx * sx;
    sx3 = sx2 * sx;

    cxp1 = (1. + cx);
    cxm1 = (1. - cx);
    cxp12 = cxp1 * cxp1;
    cxm12 = cxm1 * cxm1;
    cxm13 = cxm12 * cxm1;
    cxp13 = cxp12 * cxp1;

    temp = 0.0625 * cxp12 * cxp12;
    wigner[0] = temp;  // -4, -4 //  0

    temp = 0.0625 * cxm12 * cxm12;
    wigner[8] = temp;  //  4, -4 //  8

    temp = -0.1767766953 * cxp13 * sx;
    wigner[1] = temp;   // -3, -4 //  1
    wigner[9] = -temp;  // -4, -3 //  9

    temp = -0.1767766953 * cxm13 * sx;
    wigner[7] = temp;   //  3, -4 //  7
    wigner[17] = temp;  //  4, -3 // 17

    temp = -0.4677071733 * cxp1 * sx3;
    wigner[27] = -temp;  // -4, -1 // 27
    wigner[3] = temp;    // -1, -4 //  3

    temp = -0.4677071733 * cxm1 * sx3;
    wigner[35] = temp;  //  4, -1 // 35
    wigner[5] = temp;   //  1, -4 //  5

    temp = 0.5229125166 * sx3 * sx;
    wigner[44] = temp;  //  4,  0 // 44
    wigner[36] = temp;  // -4,  0 // 36
    wigner[4] = temp;   //  0, -4 //  4

    temp = -1.4790199458 * sx3 * cx;
    wigner[43] = temp;   //  3,  0 // 43
    wigner[37] = -temp;  // -3,  0 // 37
    wigner[13] = temp;   //  0, -3 // 13

    temp = 0.3307189139 * sx2 * cxp12;
    wigner[2] = temp;   // -2, -4 //  2
    wigner[18] = temp;  // -4, -2 // 18

    temp = 0.3307189139 * sx2 * cxm12;
    wigner[6] = temp;   //  2, -4 //  6
    wigner[26] = temp;  //  4, -2 // 26

    temp = 0.4677071733 * cxp12 * sx * (2. * cx - 1.);
    wigner[11] = -temp;  // -2, -3 // 11
    wigner[19] = temp;   // -3, -2 // 19

    temp = 0.4677071733 * cxm12 * sx * (-2. * cx - 1.);
    wigner[15] = temp;  //  2, -3 // 15
    wigner[25] = temp;  //  3, -2 // 25

    temp = 0.25 * cxp12 * (1. - 7. * cxm1 + 7. * cxm12);
    wigner[20] = temp;  // -2, -2 // 20

    temp = 0.25 * cxm12 * (1. - 7. * cxp1 + 7. * cxp12);
    wigner[24] = temp;  //  2, -2 // 24

    temp = 0.3952847075 * sx2 * (7. * cx2 - 1.0);
    wigner[42] = temp;  //  2,  0 // 42
    wigner[38] = temp;  // -2,  0 // 38
    wigner[22] = temp;  //  0, -2 // 22

    temp = 0.125 * cxp13 * (-3. + 4. * cx);
    wigner[10] = temp;  // -3, -3 // 10

    temp = 0.125 * cxm13 * (3. + 4. * cx);
    wigner[16] = temp;  //  3, -3 // 16

    temp = 0.3307189139 * cxm1 * cxp12 * (-1. + 4. * cx);
    wigner[12] = temp;  // -1, -3 // 12
    wigner[28] = temp;  // -3, -1 // 28

    temp = 0.3307189139 * cxm12 * cxp1 * (1. + 4. * cx);
    wigner[14] = temp;  //  1, -3 // 14
    wigner[34] = temp;  //  3, -1 // 34

    temp = -0.5590169944 * (4. - 18. * cxm1 + 21. * cxm12 - 7. * cxm13) * sx;
    wigner[41] = temp;   //  1,  0 // 41
    wigner[39] = -temp;  // -1,  0 // 39
    wigner[31] = temp;   //  0, -1 // 31

    temp = -0.3535533906 * (3. - 10.5 * cxm1 + 7. * cxm12) * sx * cxp1;
    wigner[29] = -temp;  // -2, -1 // 29
    wigner[21] = temp;   // -1, -2 // 21

    temp = -0.3535533906 * (10. - 17.5 * cxm1 + 7. * cxm12) * sx * cxm1;
    wigner[23] = temp;  //  1, -2 // 23
    wigner[33] = temp;  //  2, -1 // 33

    temp = 0.5 * (1. - 9. * cxm1 + 15.75 * cxm12 - 7. * cxm13) * cxp1;
    wigner[30] = temp;  // -1, -1 // 30

    temp = 0.5 * (10. - 30. * cxm1 + 26.25 * cxm12 - 7. * cxm13) * cxm1;
    wigner[32] = temp;  //  1, -1 // 32

    temp = 0.125 * (3. - 30. * cx2 + 35.0 * cx2 * cx2);
    wigner[40] = temp;  //  0,  0 // 40

    wigner += 45;  // increment counter to next matrix
  }
}

// ✅ .. note: (wigner_d_matrices_from_exp_I_beta) monitored with pytest
// ................
void wigner_d_matrices_from_exp_I_beta(const int l, const int n, const bool half,
                                       const void *restrict exp_I_beta,
                                       double *wigner) {
  if (half) {
    switch (l) {
    case 2:
      wigner_2j_half(n, exp_I_beta, wigner);
      break;
    case 4:
      wigner_4j_half(n, exp_I_beta, wigner);
      break;
    }
    return;
  }
  switch (l) {
  case 2:
    wigner_2j(n, exp_I_beta, wigner);
    break;
  case 4:
    wigner_4j(n, exp_I_beta, wigner);
    break;
  }
}

// // ✅ .. note: (__wigner_rotation) monitored with pytest
// ....................... void __wigner_rotation(const int l, const int n,
// const double *wigner,
//                        const double *cos_alpha, const complex128 *R_in,
//                        complex128 *R_out) {
//   int orientation;
//   int n1 = 2 * l + 1, m, i = 0, mp;
//   complex128 temp_initial_vector[n1], ph2, pha;
//   complex128 *final_vector;
//   int n2 = n1 * n1;

//   for (orientation = 0; orientation < n; orientation++) {

//     // calculate the alpha phase, exp(-I m alpha).
//     pha = cos_alpha[orientation] +
//           I * sqrt(1.0 - cos_alpha[orientation] * cos_alpha[orientation]);

//     // note to self: it is important to copy the pha to temporary variable,
//     ph2. ph2 = pha;

//     // copy the initial vector
//     for (m = 0; m < n1; m++) {
//       temp_initial_vector[m] = R_in[m];
//     }

//     // scale the temp initial vector with exp[-I m alpha]
//     for (m = 1; m <= l; m++) {
//       temp_initial_vector[l - m] *= ph2;
//       temp_initial_vector[l + m] *= conj(ph2);
//       ph2 *= pha;
//     }

//     final_vector = &R_out[orientation * n1];
//     i = orientation * n2;
//     // Apply wigner rotation to the temp inital vector
//     for (m = 0; m < n1; m++) {
//       final_vector[m] = 0.;
//       for (mp = 0; mp < n1; mp++) {
//         final_vector[m] += wigner[i++] * temp_initial_vector[mp];
//       }
//     }

//     // The following blas library is slower than the c for-loops
//     // probably because over head in calling blas library is more than the
//     // computation of a 5 x 5 matrix-vector multiplication.
//     // c for-loops --- 1.19 ms ± 10.4 µs
//     // blas        --- 4.19 ms ± 14.8 µs

//     // final_vector = &R_out[orientation * n1];
//     // cblas_dgemv(CblasRowMajor, CblasNoTrans, n1, n1, 1.0, &wigner[i], n1,
//     //             (double *)&temp_initial_vector[0], 2, 0.0,
//     //             (double *)&final_vector[0], 2);
//     // cblas_dgemv(CblasRowMajor, CblasNoTrans, n1, n1, 1.0, &wigner[i], n1,
//     //             (double *)&temp_initial_vector[0] + 1, 2, 0.0,
//     //             (double *)&final_vector[0] + 1, 2);
//     // i += n2;
//   }
// }

// ✅ .. note: (__wigner_rotation_2) monitored with pytest .....................
void __wigner_rotation_2(const int l, const int n, const double *wigner,
                         const void *exp_Im_alpha, const void *R_in, void *R_out) {
  double *exp_Im_alpha_ = (double *)exp_Im_alpha;
  double *R_out_ = (double *)R_out;
  double *R_in_ = (double *)R_in;

  int orientation, two_l_pm, two_l_mm;
  int n1 = 2 * l + 1, m, mp, two_l = 2 * l, two_n1 = 2 * n1;
  double S1, S2, S3, *temp, scale[5] = {1.0, -1.0, 1.0, -1.0, 1.0};
  double temp_initial_vector[18];  // maximum allocation for l=4; (2l + 1) * 2
  int location_m[5];

  for (m = 1; m <= l; m++) {
    location_m[m] = (8 - 2 * m) * n;
  }

  // Spherical tensor symmetry relations.
  // Y_{l, m} = (-1)^m Y_{l,-m}(conj)           (1)
  //
  // The following uses the above symmetry to reduce the number of calcuation steps in
  // Y_{l,m} * exp(-m I alpha) pair product for m=-l to l.
  //
  // Consider the product pairs involving +/-m, as shown below.
  //
  //    Y_{l,-m} * exp(m I alpha) = (a + ib) * (c + id)
  //                              = (ac - bd) + i(ad + bc)
  //
  // The symmetrically opposite product is then,
  //    Y_{l,m} * exp(-m I alpha) = (-1)^m (a - ib) * (c - id)
  //                              = (-1)^m (ac - bd) - i(ad + bc)
  //                              = (-1)^m (Y_{l,-m} * exp(m I alpha)) (conj)
  for (orientation = 0; orientation < n; orientation++) {
    temp_initial_vector[two_l] = R_in_[two_l];
    temp_initial_vector[two_l + 1] = R_in_[two_l + 1];

    two_l_pm = two_l;
    two_l_mm = two_l;
    for (m = 1; m <= l; m++) {
      two_l_pm += 2;
      two_l_mm -= 2;
      temp = &exp_Im_alpha_[location_m[m] + 2 * orientation];

      // temp_initial_vector[l - m] *= temp;
      // a = R_in_[two_l_mm]
      // b = R_in_[two_l_mm + 1]
      // c = *temp
      // d = *(temp + 1)
      // (a + ib) (c + id) = ac - bd + i(ad + bc)
      // Let,
      //    S1 = ac,    S2 = bd,    and   S3 = (a+b)(c+d)
      // then, (a + ib) (c + id) = S1 - S2 + i(S3 - S1 - S2)
      S1 = R_in_[two_l_mm] * *temp;
      S2 = R_in_[two_l_mm + 1] * *(temp + 1);
      S3 = (R_in_[two_l_mm] + R_in_[two_l_mm + 1]) * (*temp + *(temp + 1));
      temp_initial_vector[two_l_mm] = S1 - S2;
      temp_initial_vector[two_l_mm + 1] = S3 - S1 - S2;

      // temp_initial_vector[l + m]
      // Use symmetry relation from Eq (1).
      temp_initial_vector[two_l_pm] = temp_initial_vector[two_l_mm] * scale[m];
      temp_initial_vector[two_l_pm + 1] = -temp_initial_vector[two_l_mm + 1] * scale[m];
    }

    // Apply wigner rotation to the temp inital vector
    for (m = 0; m <= l; m++) {
      *R_out_ = 0.0;
      *(R_out_ + 1) = 0.0;
      mp = 0;
      while (mp < two_n1) {
        *R_out_ += *wigner * temp_initial_vector[mp++];
        *(R_out_ + 1) += *wigner++ * temp_initial_vector[mp++];
      }
      R_out_ += 2;
    }
  }
}

// ✅ .. note: (wigner_dm0_vector) monitored with pytest .....................
void wigner_dm0_vector(const int l, const double beta, double *R_out) {
  double sx2, sx3, cx2, cxm1, cxm12, temp;
  double cx = cos(beta), sx = sin(beta);
  switch (l) {
  case 2:
    R_out[0] = 0.6123724355 * sx * sx;  // d^2(-2,0)(beta)
    R_out[1] = 1.224744871 * sx * cx;   // d^2(-1,0)(beta)
    R_out[2] = 1.5 * cx * cx - 0.5;     // d^2(0,0)(beta)
    R_out[3] = -R_out[1];               // d^2(1,0)(beta)
    R_out[4] = R_out[0];                // d^2(2,0)(beta)
    break;
  case 4:
    sx2 = sx * sx;
    sx3 = sx2 * sx;
    cx2 = 1.0 - sx2;
    cxm1 = 1.0 - cx;
    cxm12 = cxm1 * cxm1;
    temp = 4. - 18. * cxm1 + 21. * cxm12 - 7. * cxm12 * cxm1;
    R_out[0] = 0.5229125166 * sx3 * sx;                    // d^4(-4,0)(beta)
    R_out[1] = 1.4790199458 * sx3 * cx;                    // d^4(-3,0)(beta)
    R_out[2] = 0.3952847075 * sx2 * (7. * cx2 - 1);        // d^4(-2,0)(beta)
    R_out[3] = 0.5590169944 * temp * sx;                   // d^4(-1,0)(beta)
    R_out[4] = 0.125 * (3. - 30. * cx2 + 35 * cx2 * cx2);  // d^4(0,0)(beta)
    R_out[5] = -R_out[3];                                  // d^4(1,0)(beta)
    R_out[6] = R_out[2];                                   // d^4(2,0)(beta)
    R_out[7] = -R_out[1];                                  // d^4(3,0)(beta)
    R_out[8] = R_out[0];                                   // d^4(4,0)(beta)
    break;
  }
}

/** ✅ Performs a rank l wigner rotation of the coefficients from the l rank
 * spherical tensors.
 *
 * R_out =  exp(-I m0 gamma) sum_m d^l_{m0, m} * exp(-I m alpha)
 *
 * @param l The rank of the tensor.
 * @param euler_angles A pointer to the array of three euler angles.
 * @param R_in A pointer to the array of coefficients from the l rank tensors of
 *      length 2xl+1 before rotation.
 * @param R_out A pointer to the array of coefficients from the l rank tensors
 *      of length 2xl+1 after rotation.
 */
void single_wigner_rotation(const int l, const double *euler_angles, const void *R_in,
                            void *R_out) {
  double *R_in_ = (double *)R_in;
  double *R_out_ = (double *)R_out;

  int n1 = 2 * l + 1, m, mp, k, two_l = 2 * l, two_n1 = 2 * n1;
  double real, imag, copy_real = 0.0, copy_imag = 0.0, a, b, c, d;
  double wigner[81];               // maximum allocation for l=4; (2l + 1)^2
  double temp_initial_vector[18];  // maximum allocation for l=4; (2l + 1) * 2

  // get wigner matrix corresponding to beta angle
  wigner_d_matrices(l, 1, &euler_angles[1], wigner);

  // copy the initial vector
  // cblas_zcopy(n1, R_in, 1, temp_initial_vector, 1);

  // scale the temp initial vector with exp[-I m alpha]
  // orientation at index 0, 1, 2 are alpha, beta, and gamma.
  real = cos(euler_angles[0]);
  imag = sin(euler_angles[0]);
  copy_real = real;
  copy_imag = imag;

  temp_initial_vector[two_l] = R_in_[two_l];
  temp_initial_vector[two_l + 1] = R_in_[two_l + 1];

  // scale the temp initial vector with exp[-I m alpha]
  for (m = 2; m <= two_l; m += 2) {
    // temp_initial_vector[l - m] *= exp(I m alpha); negative m
    a = R_in_[two_l - m] * real;
    b = R_in_[two_l - m + 1] * imag;
    c = R_in_[two_l - m] * imag;
    d = R_in_[two_l - m + 1] * real;
    temp_initial_vector[two_l - m] = a - b;
    temp_initial_vector[two_l - m + 1] = c + d;

    // temp_initial_vector[l + m] *= exp(I m alpha); positive m
    a = R_in_[two_l + m] * real;
    b = R_in_[two_l + m + 1] * imag;
    c = R_in_[two_l + m] * imag;
    d = R_in_[two_l + m + 1] * real;
    temp_initial_vector[two_l + m] = a + b;
    temp_initial_vector[two_l + m + 1] = -c + d;

    a = copy_real * real;
    b = copy_imag * imag;
    c = copy_real * imag;
    d = copy_imag * real;
    real = a - b;
    imag = c + d;
  }

  // for (m = 1; m <= l; m++) {
  //   temp_initial_vector[l - m] *= temp;
  //   temp_initial_vector[l + m] *= conj(temp);
  //   temp *= copy_temp;
  // }

  // Apply wigner rotation to the temp inital vector
  k = 0;
  for (m = 0; m < two_n1; m += 2) {
    R_out_[m] = 0.0;
    R_out_[m + 1] = 0.0;
    for (mp = 0; mp < two_n1; mp += 2) {
      R_out_[m] += wigner[k] * temp_initial_vector[mp];
      R_out_[m + 1] += wigner[k++] * temp_initial_vector[mp + 1];
    }
  }

  real = cos(euler_angles[2]);
  imag = sin(euler_angles[2]);
  copy_real = real;
  copy_imag = imag;
  for (m = 2; m <= two_l; m += 2) {
    // R_out_[l - m] *= exp(I m gamma); negative m
    a = R_out_[two_l - m] * real;
    b = R_out_[two_l - m + 1] * imag;
    c = R_out_[two_l - m] * imag;
    d = R_out_[two_l - m + 1] * real;
    R_out_[two_l - m] = a - b;
    R_out_[two_l - m + 1] = c + d;

    // R_out_[l + m] *= exp(-I m gamma); positive m
    a = R_out_[two_l + m] * real;
    b = R_out_[two_l + m + 1] * imag;
    c = R_out_[two_l + m] * imag;
    d = R_out_[two_l + m + 1] * real;
    R_out_[two_l + m] = a + b;
    R_out_[two_l + m + 1] = -c + d;

    a = copy_real * real;
    b = copy_imag * imag;
    c = copy_real * imag;
    d = copy_imag * real;
    real = a - b;
    imag = c + d;
  }
}

/**
 * ❌ Performs wigner rotations on a batch of wigner matrices and initial tensor
 * orientation. The wigner matrices corresponds to the beta orientations. The
 * orientations cover either an octant, hemisphere, or the sphere surface. This
 * is specified by the value of the `n_octant` variables.
 *    n_octant = 1 for octant
 *    n_octant = 4 for hemisphere
 *    n_octant = 8 for sphere
 * The function performs both second rank and fourth rank wigner rotations.
 *
 * @param octant_orientations Number of orientations on an octant.
 * @param n_octants Number of octants.
 * @param wigner_2j_matrices A pointer to a stack of 5x5 second rank wigner
 *      matrices.
 * @param R2 A pointer to the second rank tensor coefficients of length 5 to be
 *      rotated.
 * @param wigner_4j_matrices A pointer to a stack of 9x9 fourth rank wigner
 *      matrices.
 * @param R4 A pointer to the fourth rank tensor coefficients of length 9 to be
 *      rotated.
 * @param exp_Im_alpha A pointer to a `4 x octant_orientations` array with the
 *      exp(-Im α) with `octant_orientations` as the leading dimension and
 *      ordered as m=[-4,-3,-2,-1].
 * @param w2 A pointer to a stack of second rank tensor coefficients after
 *      rotation with second rank wigner matrices. The length of w2 is
 *      `octant_orientations x n_octants x 5` with 5 as the leading dimension.
 * @param w4 A pointer to a stack of fourth rank tensor coefficients after
 *      rotation with fourth rank wigner matrices. The length of w4 is
 *      `octant_orientations x n_octants x 9` with 9 as the leading dimension.
 */
void __batch_wigner_rotation(const unsigned int octant_orientations,
                             const unsigned int n_octants, double *wigner_2j_matrices,
                             complex128 *R2, double *wigner_4j_matrices, complex128 *R4,
                             complex128 *exp_Im_alpha, complex128 *w2, complex128 *w4) {
  unsigned int j, max_iter, wigner_2j_inc, wigner_4j_inc = 0;
  unsigned int w2_increment, w4_increment = 0, alpha_inc;
  double *exp_Im_alpha_ = (double *)exp_Im_alpha;

  w2_increment = 3 * octant_orientations;
  wigner_2j_inc = 5 * w2_increment;  // equal to 5 x 3 x octant_orientations;
  // 8 corresponds to m=-[-4, -3, -2, -1] times real and imag
  alpha_inc = 8 * octant_orientations;

  if (w4 != NULL) {
    w4_increment = 5 * octant_orientations;
    wigner_4j_inc = 9 * w4_increment;  // equal to 9 x 5 x octant_orientations;
  }

  max_iter = (n_octants <= 4) ? n_octants : 4;
  for (j = 0; j < max_iter; j++) {
    /* Second-rank Wigner rotation from crystal/common frame to rotor frame. */
    __wigner_rotation_2(2, octant_orientations, wigner_2j_matrices, exp_Im_alpha_, R2,
                        w2);
    w2 += w2_increment;
    if (n_octants == 8) {
      __wigner_rotation_2(2, octant_orientations, &wigner_2j_matrices[wigner_2j_inc],
                          &exp_Im_alpha_[alpha_inc], R2, w2);
      w2 += w2_increment;
    }
    if (w4 != NULL) {
      /* Fourth-rank Wigner rotation from crystal/common frame to rotor frame. */
      __wigner_rotation_2(4, octant_orientations, wigner_4j_matrices, exp_Im_alpha_, R4,
                          w4);
      w4 += w4_increment;
      if (n_octants == 8) {
        __wigner_rotation_2(4, octant_orientations, &wigner_4j_matrices[wigner_4j_inc],
                            &exp_Im_alpha_[alpha_inc], R4, w4);
        w4 += w4_increment;
      }
    }
    /**
     * Documentation
     * -------------
     *
     * Stepping the alpha phase by π/2.
     *
     * The array exp_Im_alpha is a two-dimensional array of shape
     * `4 x number_of_sidebands`, where
     *
     * exp_Im_alpha[m, :] = exp(I (4-m) alpha) for m=[0, 1, 2, 3]. When alpha += π/2
     *
     *    exp_Im_alpha[0, :] = exp(I 4alpha) * exp(I 4π/2) = exp(I 4alpha)
     *    exp_Im_alpha[0, :] *= 1
     *
     *    exp_Im_alpha[1, :] = exp(I 3alpha) * exp(I 3π/2) = -i exp(I 3alpha)
     *    exp_Im_alpha[1, :] *= -I
     *
     *    exp_Im_alpha[2, :] = exp(I 2alpha) * exp(I 2π/2) = -1 exp(I 2alpha)
     *    exp_Im_alpha[2, :] *= -1
     *
     *    exp_Im_alpha[3, :] = exp(I 1alpha) * exp(I 1π/2) = i exp(I 1alpha)
     *    exp_Im_alpha[3, :] *= I
     *
     * After four iterations, exp_Im_alpha restores to its original value.
     */
    if (n_octants != 1) {
      // m=-1
      cblas_zscal(octant_orientations, (double *)NEGATIVE_IOTA,
                  &(exp_Im_alpha_[6 * octant_orientations]), 1);
      // m=-2
      cblas_zdscal(octant_orientations, -1, &(exp_Im_alpha_[4 * octant_orientations]),
                   1);
      if (w4 != NULL) {
        // m=-3 and m=-4 rotates onto itself for octants.
        cblas_zscal(octant_orientations, (double *)IOTA,
                    &(exp_Im_alpha_[2 * octant_orientations]), 1);
      }
      if (n_octants == 8) {
        // m=-1
        cblas_zscal(octant_orientations, (double *)NEGATIVE_IOTA,
                    &(exp_Im_alpha_[14 * octant_orientations]), 1);
        // m=-2
        cblas_zdscal(octant_orientations, -1,
                     &(exp_Im_alpha_[12 * octant_orientations]), 1);
        if (w4 != NULL) {
          // m=-3 and m=-4 rotates onto itself for octants.
          cblas_zscal(octant_orientations, (double *)IOTA,
                      &(exp_Im_alpha_[10 * octant_orientations]), 1);
        }
      }
    }
  }
}

/**
 * ✅ Calculates exp(-Im angle), where angle is an array of size n.
 * The function accepts cos_angle = cos(angle)
 * The result is stored in exp_Im_angle as m x n matrix, where m = [-4,-3,-2,-1]
 * Since only negative m's are computed, the following computes exp(Im angle)
 */
void get_exp_Im_angle(const unsigned int n, const bool allow_4th_rank,
                      void *exp_Im_angle, double delta_alpha) {
  double *exp_Im_angle_ = (double *)exp_Im_angle;

  // The complex array is interpreted as alternating real and imag double array.
  // The index s_i = i*n of complex array is now at index 2*i*n.
  unsigned int m_3 = 2 * n, m_2 = 4 * n, m_1 = 6 * n;  // m_4 is 0

  // exp(-2 I angle)
  vm_double_complex_multiply(n, &exp_Im_angle_[m_1], &exp_Im_angle_[m_1],
                             &exp_Im_angle_[m_2]);

  if (allow_4th_rank) {
    // exp(-3 I angle)
    vm_double_complex_multiply(n, &exp_Im_angle_[m_2], &exp_Im_angle_[m_1],
                               &exp_Im_angle_[m_3]);
    // exp(-4 I angle)
    vm_double_complex_multiply(n, &exp_Im_angle_[m_3], &exp_Im_angle_[m_1],
                               exp_Im_angle_);
  }

  if (delta_alpha != 0.0) {
    double *exp_da = malloc_double(8);
    unsigned int m_4p = 8 * n, m_3p = 10 * n, m_2p = 12 * n, m_1p = 14 * n;

    exp_da[0] = cos(delta_alpha);
    exp_da[1] = sin(delta_alpha);
    exp_da[2] = cos(2.0 * delta_alpha);
    exp_da[3] = sin(2.0 * delta_alpha);
    exp_da[4] = cos(3.0 * delta_alpha);
    exp_da[5] = sin(3.0 * delta_alpha);
    exp_da[6] = cos(4.0 * delta_alpha);
    exp_da[7] = sin(4.0 * delta_alpha);

    // exp(-I angle) * exp(-I delta angle)
    cblas_zcopy(n, &exp_Im_angle_[m_1], 1, &exp_Im_angle_[m_1p], 1);
    cblas_zscal(n, exp_da, &exp_Im_angle_[m_1p], 1);

    // exp(-2I angle) * exp(-2I delta angle)
    cblas_zcopy(n, &exp_Im_angle_[m_2], 1, &exp_Im_angle_[m_2p], 1);
    cblas_zscal(n, &exp_da[2], &exp_Im_angle_[m_2p], 1);

    if (allow_4th_rank) {
      // exp(-3I angle) * exp(-3I delta angle)
      cblas_zcopy(n, &exp_Im_angle_[m_3], 1, &exp_Im_angle_[m_3p], 1);
      cblas_zscal(n, &exp_da[4], &exp_Im_angle_[m_3p], 1);

      // exp(-4I angle) * exp(-4I delta angle)
      cblas_zcopy(n, exp_Im_angle_, 1, &exp_Im_angle_[m_4p], 1);
      cblas_zscal(n, &exp_da[6], &exp_Im_angle_[m_4p], 1);
    }
  }
}
