//
//  angular_momentum.c
//
//  Created by Philip Grandinetti on 4/12/17.
//  Copyright © 2017 Philip Grandinetti. All rights reserved.
//  Contribution: Deepansh J. Srivatava. contact: srivastava.89@osu.edu
//      Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com

#include "angular_momentum.h"

/* calculate Wigner rotation matrices */

/* This routine calculates the factorial of x */

double fac(double x) {
  double sum = 1;
  int ix;

  if (x < 0) {
    exit(1);
  }
  ix = (int)x;
  for (; ix > 1; ix--) {
    sum *= ix;
  }
  return sum;
}

/* power function */
double my_power(double x, int n) {
  double temp;
  if (n == 0) {
    return (1.);
  }
  temp = 1.;
  for (; n >= 1; n--) {
    temp *= x;
  }
  return (temp);
}

// ✅ .. note: (wigner_d_matrix) monitored with pytest .........................
void wigner_d_matrix(const int l, const int n, const double *angle,
                     double *wigner) {
  double *cos_angle = malloc_double(n);
  vm_double_cosine(n, angle, cos_angle);
  __wigner_d_matrix_cosine(l, n, cos_angle, wigner);
  free(cos_angle);
}

// ✅ .. note: (__wigner_d_matrix_cosine) monitored with pytest ................
void __wigner_d_matrix_cosine(const int l, const int n,
                              const double *restrict cos_angle,
                              double *wigner) {
  double cx, cx2, sx, temp;
  int i;
  if (l == 2) {
    double t1;
    for (i = 0; i < n; i++) {
      cx = *cos_angle++;
      cx2 = cx * cx;
      sx = sqrt(1. - cx2);

      t1 = (1. + cx);
      temp = -sx * t1 * 0.5;
      wigner[19] = temp;  //  2,  1 // 19
      wigner[5] = -temp;  // -2, -1 //  5
      wigner[23] = -temp; //  1,  2 // 23
      wigner[1] = temp;   // -1, -2 //  1

      temp = t1 * t1 * 0.25;
      wigner[24] = temp; //  2,  2 // 24
      wigner[0] = temp;  // -2, -2 //  0

      t1 = (1. - cx);
      temp = -sx * t1 * 0.5;
      wigner[9] = temp;   //  2, -1 //  9
      wigner[15] = -temp; // -2,  1 // 15
      wigner[3] = temp;   //  1, -2 //  3
      wigner[21] = -temp; // -1,  2 // 21

      temp = t1 * t1 * 0.25;
      wigner[4] = temp;  //  2, -2 //  4
      wigner[20] = temp; // -2,  2 // 20

      temp = 0.6123724355 * sx * sx;
      wigner[14] = temp; //  2,  0 // 14
      wigner[10] = temp; // -2,  0 // 10
      wigner[22] = temp; //  0,  2 // 22
      wigner[2] = temp;  //  0, -2 //  2

      temp = 1.224744871 * sx * cx;
      wigner[13] = -temp; //  1,  0 // 13
      wigner[17] = temp;  //  0,  1 // 17
      wigner[7] = -temp;  //  0, -1 //  7
      wigner[11] = temp;  // -1,  0 // 11

      temp = (2.0 * cx2 + cx - 1.) * 0.5;
      wigner[18] = temp; //  1,  1 // 18
      wigner[6] = temp;  // -1, -1 //  6

      temp = -(2.0 * cx2 - cx - 1.) * 0.5;
      wigner[8] = temp;  //  1, -1 //  8
      wigner[16] = temp; // -1,  1 // 16

      wigner[12] = 1.5 * cx2 - 0.5; // 0,  0 // 12

      wigner += 25;
    }
  }
  if (l == 4) {
    double sx2, sx3, cxp1, cxm1, cxp12, cxm12, cxm13, cxp13;
    for (i = 0; i < n; i++) {
      cx = *cos_angle++;
      cx2 = cx * cx;
      sx = sqrt(1. - cx2);

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
      wigner[80] = temp; //  4,  4 // 80

      temp = 0.0625 * cxm12 * cxm12;
      wigner[72] = temp; // -4,  4 // 72
      wigner[8] = temp;  //  4, -4 //  8

      temp = -0.1767766953 * cxp13 * sx;
      wigner[1] = temp;   // -3, -4 //  1
      wigner[9] = -temp;  // -4, -3 //  9
      wigner[79] = -temp; //  3,  4 // 79
      wigner[71] = temp;  //  4,  3 //  9

      temp = -0.1767766953 * cxm13 * sx;
      wigner[7] = temp;   //  3, -4 //  7
      wigner[63] = -temp; // -4,  3 // 63
      wigner[73] = -temp; // -3,  4 // 73
      wigner[17] = temp;  //  4, -3 // 17

      temp = -0.4677071733 * cxp1 * sx3;
      wigner[53] = temp;  //  4,  1 // 53
      wigner[27] = -temp; // -4, -1 // 27
      wigner[77] = -temp; //  1,  4 // 77
      wigner[3] = temp;   // -1, -4 //  3

      temp = -0.4677071733 * cxm1 * sx3;
      wigner[35] = temp;  //  4, -1 // 35
      wigner[45] = -temp; // -4,  1 // 45
      wigner[75] = -temp; // -1,  4 // 75
      wigner[5] = temp;   //  1, -4 //  5

      temp = 0.5229125166 * sx3 * sx;
      wigner[44] = temp; //  4,  0 // 44
      wigner[36] = temp; // -4,  0 // 36
      wigner[76] = temp; //  0,  4 // 76
      wigner[4] = temp;  //  0, -4 //  4

      temp = -1.4790199458 * sx3 * cx;
      wigner[43] = temp;  //  3,  0 // 43
      wigner[37] = -temp; // -3,  0 // 37
      wigner[67] = -temp; //  0,  3 // 67
      wigner[13] = temp;  //  0, -3 // 13

      temp = 0.3307189139 * sx2 * cxp12;
      wigner[78] = temp; //  2,  4 // 78
      wigner[2] = temp;  // -2, -4 //  2
      wigner[62] = temp; //  4,  2 // 62
      wigner[18] = temp; // -4, -2 // 18

      temp = 0.3307189139 * sx2 * cxm12;
      wigner[6] = temp;  //  2, -4 //  6
      wigner[74] = temp; // -2,  4 // 74
      wigner[54] = temp; // -4,  2 // 54
      wigner[26] = temp; //  4, -2 // 26

      temp = 0.4677071733 * cxp12 * sx * (2. * cx - 1.);
      wigner[69] = temp;  //  2,  3 // 69
      wigner[11] = -temp; // -2, -3 // 11
      wigner[61] = -temp; //  3,  2 // 61
      wigner[19] = temp;  // -3, -2 // 19

      temp = 0.4677071733 * cxm12 * sx * (-2. * cx - 1.);
      wigner[15] = temp;  //  2, -3 // 15
      wigner[65] = -temp; // -2,  3 // 65
      wigner[55] = -temp; // -3,  2 // 55
      wigner[25] = temp;  //  3, -2 // 25

      temp = 0.25 * cxp12 * (1. - 7. * cxm1 + 7. * cxm12);
      wigner[60] = temp; //  2,  2 // 60
      wigner[20] = temp; // -2, -2 // 20

      temp = 0.25 * cxm12 * (1. - 7. * cxp1 + 7. * cxp12);
      wigner[56] = temp; // -2,  2 // 56
      wigner[24] = temp; //  2, -2 // 24

      temp = 0.3952847075 * sx2 * (7. * cx2 - 1.0);
      wigner[42] = temp; //  2,  0 // 42
      wigner[38] = temp; // -2,  0 // 38
      wigner[58] = temp; //  0,  2 // 58
      wigner[22] = temp; //  0, -2 // 22

      temp = 0.125 * cxp13 * (-3. + 4. * cx);
      wigner[10] = temp; // -3, -3 // 10
      wigner[70] = temp; //  3,  3 // 70

      temp = 0.125 * cxm13 * (3. + 4. * cx);
      wigner[64] = temp; // -3,  3 // 64
      wigner[16] = temp; //  3, -3 // 16

      temp = 0.3307189139 * cxm1 * cxp12 * (-1. + 4. * cx);
      wigner[12] = temp; // -1, -3 // 12
      wigner[28] = temp; // -3, -1 // 28
      wigner[68] = temp; //  1,  3 // 68
      wigner[52] = temp; //  3,  1 // 52

      temp = 0.3307189139 * cxm12 * cxp1 * (1. + 4. * cx);
      wigner[14] = temp; //  1, -3 // 14
      wigner[46] = temp; // -3,  1 // 46
      wigner[66] = temp; // -1,  3 // 66
      wigner[34] = temp; //  3, -1 // 34

      temp = -0.5590169944 * (4. - 18. * cxm1 + 21. * cxm12 - 7. * cxm13) * sx;
      wigner[41] = temp;  //  1,  0 // 41
      wigner[39] = -temp; // -1,  0 // 39
      wigner[49] = -temp; //  0,  1 // 49
      wigner[31] = temp;  //  0, -1 // 31

      temp = -0.3535533906 * (3. - 10.5 * cxm1 + 7. * cxm12) * sx * cxp1;
      wigner[51] = temp;  //  2,  1 // 51
      wigner[29] = -temp; // -2, -1 // 29
      wigner[59] = -temp; //  1,  2 // 59
      wigner[21] = temp;  // -1, -2 // 21

      temp = -0.3535533906 * (10. - 17.5 * cxm1 + 7. * cxm12) * sx * cxm1;
      wigner[23] = temp;  //  1, -2 // 23
      wigner[57] = -temp; // -1,  2 // 57
      wigner[47] = -temp; // -2,  1 // 47
      wigner[33] = temp;  //  2, -1 // 33

      temp = 0.5 * (1. - 9. * cxm1 + 15.75 * cxm12 - 7. * cxm13) * cxp1;
      wigner[30] = temp; // -1, -1 // 30
      wigner[50] = temp; //  1,  1 // 50

      temp = 0.5 * (10. - 30. * cxm1 + 26.25 * cxm12 - 7. * cxm13) * cxm1;
      wigner[32] = temp; //  1, -1 // 32
      wigner[48] = temp; // -1,  1 // 48

      temp = 0.125 * (3. - 30. * cx2 + 35.0 * cx2 * cx2);
      wigner[40] = temp; //  0,  0 // 40

      wigner += 81;
    }
  }
}

// ✅ .. note: (__wigner_rotation) monitored with pytest .......................
void __wigner_rotation(const int l, const int n, const double *wigner,
                       const double *cos_alpha, const complex128 *R_in,
                       complex128 *R_out) {
  int orientation;
  int n1 = 2 * l + 1, m, i = 0, mp;
  complex128 temp_initial_vector[n1], ph2, pha;
  complex128 *final_vector;
  int n2 = n1 * n1;

  for (orientation = 0; orientation < n; orientation++) {

    // calculate the alpha phase, exp(-I m alpha).
    pha = cos_alpha[orientation] +
          I * sqrt(1.0 - cos_alpha[orientation] * cos_alpha[orientation]);

    // note to self: it is important to copy the pha to temporary variable, ph2.
    ph2 = pha;

    // copy the initial vector
    for (m = 0; m < n1; m++) {
      temp_initial_vector[m] = R_in[m];
    }

    // scale the temp initial vector with exp[-I m alpha]
    for (m = 1; m <= l; m++) {
      temp_initial_vector[l - m] *= ph2;
      temp_initial_vector[l + m] *= conj(ph2);
      ph2 *= pha;
    }

    final_vector = &R_out[orientation * n1];
    i = orientation * n2;
    // Apply wigner rotation to the temp inital vector
    for (m = 0; m < n1; m++) {
      final_vector[m] = 0.;
      for (mp = 0; mp < n1; mp++) {
        final_vector[m] += wigner[i++] * temp_initial_vector[mp];
      }
    }

    // The following blas library is slower than the c for-loops
    // probably because over head in calling blas library is more than the
    // computation of a 5 x 5 matrix-vector multiplication.
    // c for-loops --- 1.19 ms ± 10.4 µs
    // blas        --- 4.19 ms ± 14.8 µs

    // final_vector = &R_out[orientation * n1];
    // cblas_dgemv(CblasRowMajor, CblasNoTrans, n1, n1, 1.0, &wigner[i], n1,
    //             (double *)&temp_initial_vector[0], 2, 0.0,
    //             (double *)&final_vector[0], 2);
    // cblas_dgemv(CblasRowMajor, CblasNoTrans, n1, n1, 1.0, &wigner[i], n1,
    //             (double *)&temp_initial_vector[0] + 1, 2, 0.0,
    //             (double *)&final_vector[0] + 1, 2);
    // i += n2;
  }
}

#if !(__STDC_VERSION__ >= 199901L)
void __wigner_rotation_2(const int l, const int n, const double *wigner,
                         const complex128 *exp_Im_alpha, const complex128 *R_in,
                         complex128 *R_out) {
  int orientation;
  int n1 = 2 * l + 1, m, i = 0, mp;
  complex128 *temp_initial_vector = malloc_complex128(n1);
  complex128 *final_vector, temp;

  for (orientation = 0; orientation < n; orientation++) {

    // copy the initial vector
    for (m = 0; m < n1; m++) {
      temp_initial_vector[m] = R_in[m];
    }

    // scale the temp initial vector with exp[-I m alpha]
    for (m = 1; m <= l; m++) {
      temp = exp_Im_alpha[(4 - m) * n + orientation];
      temp_initial_vector[l - m] =
          complex_multiply(temp_initial_vector[l - m], temp);
      temp_initial_vector[l + m] =
          complex_multiply(temp_initial_vector[l + m], conj(temp));
    }

    final_vector = &R_out[orientation * n1];
    // Apply wigner rotation to the temp inital vector
    for (m = 0; m < n1; m++) {
      final_vector[m] = 0.;
      for (mp = 0; mp < n1; mp++) {
        *final_vector.real += *wigner * temp_initial_vector[mp].real;
        *final_vector.imag += *wigner++ * temp_initial_vector[mp].imag;
      }
      final_vector++;
    }
  }
}
#else

// ✅ .. note: (__wigner_rotation_2) monitored with pytest .....................
void __wigner_rotation_2(const int l, const int n, const double *wigner,
                         const complex128 *exp_Im_alpha, const complex128 *R_in,
                         complex128 *R_out) {
  int orientation;
  int n1 = 2 * l + 1, m, mp;
  complex128 *temp_initial_vector = malloc_complex128(n1);
  complex128 temp;

  for (orientation = 0; orientation < n; orientation++) {

    // copy the initial vector
    cblas_zcopy(n1, R_in, 1, temp_initial_vector, 1);

    // scale the temp initial vector with exp[-I m alpha]
    for (m = 1; m <= l; m++) {
      temp = exp_Im_alpha[(4 - m) * n + orientation];
      temp_initial_vector[l - m] *= temp;
      temp_initial_vector[l + m] *= conj(temp);
    }

    // Apply wigner rotation to the temp inital vector
    for (m = 0; m < n1; m++) {
      *R_out = 0.0;
      for (mp = 0; mp < n1; mp++) {
        *R_out += *wigner++ * temp_initial_vector[mp];
      }
      R_out++;
    }
  }
}
#endif

// ✅ .. note: (__wigner_dm0_vector) monitored with pytest .....................
void __wigner_dm0_vector(const int l, const double beta, double *R_out) {
  double cx = cos(beta), sx = sin(beta);
  if (l == 2) {
    R_out[0] = 0.6123724355 * sx * sx;
    R_out[1] = 1.224744871 * sx * cx;
    R_out[2] = 1.5 * cx * cx - 0.5;
    R_out[3] = -R_out[1];
    R_out[4] = R_out[0];
  }
  if (l == 4) {
    double sx2 = sx * sx, sx3 = sx2 * sx, cx2 = 1.0 - sx2;
    double cxm1 = 1.0 - cx, cxm12 = cxm1 * cxm1;
    double temp = 4. - 18. * cxm1 + 21. * cxm12 - 7. * cxm12 * cxm1;
    R_out[0] = 0.5229125166 * sx3 * sx;
    R_out[1] = 1.4790199458 * sx3 * cx;
    R_out[2] = 0.3952847075 * sx2 * (7. * cx2 - 1);
    R_out[3] = 0.5590169944 * temp * sx;
    R_out[4] = 0.125 * (3. - 30. * cx2 + 35 * cx2 * cx2);
    R_out[5] = -R_out[3];
    R_out[6] = R_out[2];
    R_out[7] = -R_out[1];
    R_out[8] = R_out[0];
  }
}

// ❌.. note: (full_DLM) not tested ............................................
void full_DLM(complex128 *wigner, int l, double *omega) {
  if (l == 2) {
    complex128 pha[5], phg[5];
    int m;
    for (m = -2; m <= 2; m++) {
      pha[m + 2] = cexp(-I * m * omega[0]);
      // printf("pha %f \n", creal(pha[m+2]));
      phg[m + 2] = cexp(-I * m * omega[2]);
    }
    double cx = cos(omega[1]);
    double sx = sin(omega[1]);

    double t1 = (1. + cx);
    double temp = -sx * t1 / 2.;
    wigner[19] = temp * pha[4] * phg[3];  //  2,  1  // 19
    wigner[5] = -temp * pha[0] * phg[1];  // -2, -1 // 5
    wigner[23] = -temp * pha[3] * phg[4]; //  1,  2 // 23
    wigner[1] = temp * pha[1] * phg[0];   // -1, -2 // 1

    temp = t1 * t1 / 4.;
    wigner[24] = temp * pha[4] * phg[4]; //  2,  2 // 24
    wigner[0] = temp * pha[0] * phg[0];  // -2, -2 // 0

    t1 = (1. - cx);
    temp = -sx * t1 / 2.;
    wigner[9] = temp * pha[4] * phg[1];   //  2, -1 // 9
    wigner[15] = -temp * pha[0] * phg[3]; // -2,  1 // 15
    wigner[3] = temp * pha[3] * phg[0];   //  1, -2 // 3
    wigner[21] = -temp * pha[1] * phg[4]; // -1,  2 // 21

    temp = t1 * t1 / 4.;
    wigner[4] = temp * pha[4] * phg[0];  //  2, -2 // 4
    wigner[20] = temp * pha[0] * phg[4]; // -2,  2 // 20

    temp = 0.6123724355 * sx * sx;
    wigner[14] = temp * pha[4] * phg[2]; //  2,  0 //14
    wigner[10] = temp * pha[0] * phg[2]; // -2,  0 // 10
    wigner[22] = temp * pha[2] * phg[4]; //  0,  2 // 22
    wigner[2] = temp * pha[2] * phg[0];  //  0, -2 // 2

    temp = 1.224744871 * sx * cx;
    wigner[13] = -temp * pha[3] * phg[2]; //  1,  0 // 13
    wigner[17] = temp * pha[2] * phg[3];  //  0,  1 // 17
    wigner[7] = -temp * pha[2] * phg[1];  //  0, -1 // 7
    wigner[11] = temp * pha[1] * phg[2];  // -1,  0 // 11

    temp = (2.0 * cx * cx + cx - 1.) / 2.;
    wigner[18] = temp * pha[3] * phg[3]; //  1,  1 //18
    wigner[6] = temp * pha[1] * phg[1];  // -1, -1 // 6

    temp = -(2.0 * cx * cx - cx - 1.) / 2.;
    wigner[8] = temp * pha[3] * phg[1];  //  1, -1 // 8
    wigner[16] = temp * pha[1] * phg[3]; // -1,  1 // 16

    wigner[12] = 1.5 * cx * cx - .5 * pha[2] * phg[2]; // 0,  0 // 12
  }
  if (l == 4) {
    complex128 pha[9]; //, phg[9];
    int m;
    for (m = -4; m <= 4; m++) {
      pha[m + 4] = cexp(-I * m * omega[0]);
      // printf("pha %f \n", creal(pha[m+2]));
      // phg[m + 4] = cexp(-I * m * omega[2]);
    }
    double cx = cos(omega[1]);
    double sx = sin(omega[1]);
    double sx2 = sx * sx, sx3 = sx2 * sx;
    double cx2 = cx * cx;

    double cxp1 = (1. + cx), cxm1 = (1. - cx), cxp12 = cxp1 * cxp1,
           cxm12 = cxm1 * cxm1;

    // index = (gamma+4)*9 + (alpha+4)
    double temp = 0.0625 * cxp12 * cxp12;
    wigner[0] = temp * pha[0];  // -4, -4 //  0
    wigner[80] = temp * pha[8]; //  4,  4 // 80

    temp = 0.0625 * cxm12 * cxm12;
    wigner[72] = temp * pha[0]; // -4,  4 // 72
    wigner[8] = temp * pha[8];  //  4, -4 //  8

    temp = -0.1767766953 * cxp12 * cxp1 * sx;
    wigner[1] = temp * pha[1];   // -3, -4 //  1
    wigner[9] = -temp * pha[0];  // -4, -3 //  9
    wigner[79] = -temp * pha[7]; //  3,  4 // 79
    wigner[71] = temp * pha[8];  //  4,  3 // 71

    temp = -0.1767766953 * cxm12 * cxm1 * sx;
    wigner[7] = temp * pha[7];   //  3, -4 //  7
    wigner[63] = -temp * pha[0]; // -4,  3 // 63
    wigner[73] = -temp * pha[1]; // -3,  4 // 73
    wigner[17] = temp * pha[8];  //  4, -3 // 17

    temp = -0.4677071733 * cxp1 * sx3;
    wigner[53] = temp * pha[8];  //  4,  1 // 53
    wigner[27] = -temp * pha[0]; // -4, -1 // 27
    wigner[77] = -temp * pha[5]; //  1,  4 // 77
    wigner[3] = temp * pha[3];   // -1, -4 //  3

    temp = -0.4677071733 * cxm1 * sx3;
    wigner[35] = temp * pha[8];  //  4, -1 // 35
    wigner[45] = -temp * pha[0]; // -4,  1 // 45
    wigner[75] = -temp * pha[3]; // -1,  4 // 75
    wigner[5] = temp * pha[5];   //  1, -4 //  5

    temp = 0.5229125166 * sx3 * sx;
    wigner[44] = temp * pha[8]; //  4,  0 // 44
    wigner[36] = temp * pha[0]; // -4,  0 // 36
    wigner[76] = temp * pha[4]; //  0,  4 // 76
    wigner[4] = temp * pha[4];  //  0, -4 //  4

    temp = -1.4790199458 * sx3 * cx;
    wigner[43] = temp * pha[7];  //  3,  0 // 43
    wigner[37] = -temp * pha[1]; // -3,  0 // 37
    wigner[67] = -temp * pha[4]; //  0,  3 // 67
    wigner[13] = temp * pha[4];  //  0, -3 // 13

    temp = 0.3307189139 * sx2 * cxp12;
    wigner[78] = temp * pha[6]; //  2,  4 // 78
    wigner[2] = temp * pha[2];  // -2, -4 //  2
    wigner[62] = temp * pha[8]; //  4,  2 // 62
    wigner[18] = temp * pha[0]; // -4, -2 // 18

    temp = 0.3307189139 * sx2 * cxm12;
    wigner[6] = temp * pha[6];  //  2, -4 //  6
    wigner[74] = temp * pha[2]; // -2,  4 // 74
    wigner[54] = temp * pha[0]; // -4,  2 // 54
    wigner[26] = temp * pha[8]; //  4, -2 // 26

    temp = 0.4677071733 * cxp12 * sx * (2. * cx - 1.);
    wigner[69] = temp * pha[6];  //  2,  3 // 69
    wigner[11] = -temp * pha[2]; // -2, -3 // 11
    wigner[61] = -temp * pha[7]; //  3,  2 // 61
    wigner[19] = temp * pha[1];  // -3, -2 // 19

    temp = 0.4677071733 * cxm12 * sx * (-2. * cx - 1.);
    wigner[15] = temp * pha[6];  //  2, -3 // 15
    wigner[65] = -temp * pha[2]; // -2,  3 // 65
    wigner[55] = -temp * pha[1]; // -3,  2 // 55
    wigner[25] = temp * pha[7];  //  3, -2 // 25

    temp = 0.25 * cxp12 * (1. - 7. * cxm1 + 7. * cxm12);
    wigner[60] = temp * pha[6]; //  2,  2 // 60
    wigner[20] = temp * pha[2]; // -2, -2 // 20

    temp = 0.25 * cxm12 * (1. - 7. * cxp1 + 7. * cxp12);
    wigner[56] = temp * pha[2]; // -2,  2 // 56
    wigner[24] = temp * pha[6]; //  2, -2 // 24

    temp = 0.3952847075 * sx2 * (7. * cx2 - 1);
    wigner[42] = temp * pha[6]; //  2,  0 // 42
    wigner[38] = temp * pha[2]; // -2,  0 // 38
    wigner[58] = temp * pha[4]; //  0,  2 // 58
    wigner[22] = temp * pha[4]; //  0, -2 // 22

    temp = 0.125 * cxp12 * cxp1 * (-3. + 4. * cx);
    wigner[10] = temp * pha[1]; // -3, -3 // 10
    wigner[70] = temp * pha[7]; //  3,  3 // 70

    temp = 0.125 * cxm12 * cxm1 * (3. + 4. * cx);
    wigner[64] = temp * pha[1]; // -3,  3 // 64
    wigner[16] = temp * pha[7]; //  3, -3 // 16

    temp = 0.3307189139 * cxm1 * cxp12 * (-1. + 4. * cx);
    wigner[12] = temp * pha[3]; // -1, -3 // 12
    wigner[28] = temp * pha[1]; // -3, -1 // 28
    wigner[68] = temp * pha[5]; //  1,  3 // 68
    wigner[52] = temp * pha[7]; //  3,  1 // 52

    temp = 0.3307189139 * cxm12 * cxp1 * (1. + 4. * cx);
    wigner[14] = temp * pha[5]; //  1, -3 // 14
    wigner[46] = temp * pha[1]; // -3,  1 // 46
    wigner[66] = temp * pha[3]; // -1,  3 // 66
    wigner[34] = temp * pha[7]; //  3, -1 // 34

    temp = -0.5590169944 * (4. - 18. * cxm1 + 21. * cxm12 - 7. * cxm12 * cxm1) *
           sx;
    wigner[41] = temp * pha[5];  //  1,  0 // 41
    wigner[39] = -temp * pha[3]; // -1,  0 // 39
    wigner[49] = -temp * pha[4]; //  0,  1 // 49
    wigner[31] = temp * pha[4];  //  0, -1 // 31

    temp = -0.3535533906 * (3. - 10.5 * cxm1 + 7. * cxm12) * sx * cxp1;
    wigner[51] = temp * pha[6];  //  2,  1 // 51
    wigner[29] = -temp * pha[2]; // -2, -1 // 29
    wigner[59] = -temp * pha[5]; //  1,  2 // 59
    wigner[21] = temp * pha[3];  // -1, -2 // 21

    temp = -0.3535533906 * (10. - 17.5 * cxm1 + 7. * cxm12) * sx * cxm1;
    wigner[23] = temp * pha[5];  //  1, -2 // 23
    wigner[57] = -temp * pha[3]; // -1,  2 // 57
    wigner[47] = -temp * pha[2]; // -2,  1 // 47
    wigner[33] = temp * pha[6];  //  2, -1 // 33

    temp = 0.5 * (1. - 9. * cxm1 + 15.75 * cxm12 - 7. * cxm12 * cxm1) * cxp1;
    wigner[30] = temp * pha[3]; // -1, -1 // 30
    wigner[50] = temp * pha[5]; //  1,  1 // 50

    temp = 0.5 * (10. - 30. * cxm1 + 26.25 * cxm12 - 7. * cxm12 * cxm1) * cxm1;
    wigner[32] = temp * pha[5]; //  1, -1 // 32
    wigner[48] = temp * pha[3]; // -1,  1 // 48

    temp = 0.125 * (3. - 30. * cx2 + 35. * cx2 * cx2);
    wigner[40] = temp * pha[4]; //  0,  0 // 40
  }
}

// ❌.. note: (full_DLM_trig) not tested .......................................
// @parameters
// int l: The angular momentum quantum number
// complex128 *wigner: a pointer to (2l+1)*(2l+1) size matrix of full
// Wigner matrix
void full_DLM_trig(complex128 *wigner, int l, double cosAlpha, double sinAlpha,
                   double cosBeta, double sinBeta) {
  //    double cosGamma,
  //    double sinGamma){
  if (l == 2) {
    // pha[m+2] holds the value of exp(-I m alpha)
    // phg[m+2] holds the value of exp(-I m beta)
    complex128 pha[5]; //, phg[5];
    pha[2] = 1.0;
    // phg[2] = 1.0;

    // for m=-1
    // ph[1] = cos(-1 alpha) - I sin(-1 alpha)
    //       = cos(alpha) + I sin(alpha)
    pha[1] = cosAlpha + I * sinAlpha;
    // phg[1] = cosGamma + I*sinGamma;

    // similarly,
    // for m=1
    // ph[3] = cos(1 alpha) - I sin(1 alpha)
    //       = cos(alpha) - I sin(alpha)
    pha[3] = cosAlpha - I * sinAlpha;
    // phg[3] = cosGamma - I*sinGamma;

    // for m= -2
    // pha[0] = cos(-2 alpha) - I sin(-2 alpha)
    //        = cos(alpha)^2 - sin(alpha)^2 + I 2 sin(alpha)cos(alpha)
    //        = ( cos(alpha) + I sin(alpha) )^2
    //        = (ph[1])^2
    pha[0] = pha[1] * pha[1];
    // phg[0] = cpow(phg[1], 2);

    // similarly for m=2
    // pha[4] = cos(2 alpha) - I sin(2 alpha)
    //        = cos(alpha)^2 - sin(alpha)^2 - I 2 sin(alpha)cos(alpha)
    //        = ( cos(alpha) - I sin(alpha) )^2
    //        = (ph[3])^2
    pha[4] = pha[3] * pha[3];
    // phg[4] = cpow(phg[3], 2);

    double cx = cosBeta;
    double sx = sinBeta;
    double cx2 = cx * cx;
    // calculating wigner small d terms.

    double t1 = (1. + cx);
    double temp = -sx * t1 / 2.;
    wigner[19] = temp * pha[4];  // * phg[3];    //  2,  1 // 19
    wigner[5] = -temp * pha[0];  // * phg[1];    // -2, -1 //  5
    wigner[23] = -temp * pha[3]; // * phg[4];   //  1,  2 // 23
    wigner[1] = temp * pha[1];   // * phg[0];     // -1, -2 //  1

    temp = t1 * t1 / 4.;
    wigner[24] = temp * pha[4]; // * phg[4];    //  2,  2 // 24
    wigner[0] = temp * pha[0];  // * phg[0];     // -2, -2 //  0

    t1 = (1. - cx);
    temp = -sx * t1 / 2.;
    wigner[9] = temp * pha[4];   // * phg[1];     //  2, -1 //  9
    wigner[15] = -temp * pha[0]; // * phg[3];   // -2,  1 // 15
    wigner[3] = temp * pha[3];   // * phg[0];     //  1, -2 //  3
    wigner[21] = -temp * pha[1]; // * phg[4];   // -1,  2 // 21

    temp = t1 * t1 / 4.;
    wigner[4] = temp * pha[4];  // * phg[0];     //  2, -2 //  4
    wigner[20] = temp * pha[0]; // * phg[4];    // -2,  2 // 20

    temp = 0.6123724355 * sx * sx;
    wigner[14] = temp * pha[4]; // * phg[2];    //  2,  0 // 14
    wigner[10] = temp * pha[0]; // * phg[2];    // -2,  0 // 10
    wigner[22] = temp * pha[2]; // * phg[4];    //  0,  2 // 22
    wigner[2] = wigner[22]; // temp * pha[2];// * phg[0];     //  0, -2 //  2

    temp = 1.224744871 * sx * cx;
    wigner[13] = -temp * pha[3]; // * phg[2];   //  1,  0 // 13
    wigner[17] = temp * pha[2];  // * phg[3];    //  0,  1 // 17
    wigner[7] = -wigner[17];    //-temp * pha[2];// * phg[1];    //  0, -1 //  7
    wigner[11] = temp * pha[1]; // * phg[2];    // -1,  0 // 11

    temp = (2.0 * cx2 + cx - 1.) / 2.;
    wigner[18] = temp * pha[3]; // * phg[3];    //  1,  1 // 18
    wigner[6] = temp * pha[1];  // * phg[1];     // -1, -1 //  6

    temp = -(2.0 * cx2 - cx - 1.) / 2.;
    wigner[8] = temp * pha[3];  // * phg[1];     //  1, -1 //  8
    wigner[16] = temp * pha[1]; // * phg[3];    // -1,  1 // 16

    // 0,  0 // 12
    wigner[12] = 1.5 * cx2 - .5 * pha[2]; // * phg[2];
  }

  if (l == 4) {
    complex128 pha[9]; //, phg[5];
    pha[4] = 1.0;
    pha[5] = cosAlpha - I * sinAlpha;
    pha[3] = cosAlpha + I * sinAlpha;
    pha[6] = pha[5] * pha[5];
    pha[2] = conj(pha[6]);
    pha[7] = pha[6] * pha[5];
    pha[1] = conj(pha[7]);
    pha[8] = pha[7] * pha[5];
    pha[0] = conj(pha[8]);

    double cx = cosBeta;
    double sx = sinBeta;
    double sx2 = sx * sx, sx3 = sx2 * sx;
    double cx2 = cx * cx;

    double cxp1 = (1. + cx), cxm1 = (1. - cx), cxp12 = cxp1 * cxp1,
           cxm12 = cxm1 * cxm1;

    // index = (gamma+4)*9 + (alpha+4)
    double temp = 0.0625 * cxp12 * cxp12;
    wigner[0] = temp * pha[0];  // -4, -4 //  0
    wigner[80] = temp * pha[8]; //  4,  4 // 80

    temp = 0.0625 * cxm12 * cxm12;
    wigner[72] = temp * pha[0]; // -4,  4 // 72
    wigner[8] = temp * pha[8];  //  4, -4 //  8

    temp = -0.1767766953 * cxp12 * cxp1 * sx;
    wigner[1] = temp * pha[1];   // -3, -4 //  1
    wigner[9] = -temp * pha[0];  // -4, -3 //  9
    wigner[79] = -temp * pha[7]; //  3,  4 // 79
    wigner[71] = temp * pha[8];  //  4,  3 // 71

    temp = -0.1767766953 * cxm12 * cxm1 * sx;
    wigner[7] = temp * pha[7];   //  3, -4 //  7
    wigner[63] = -temp * pha[0]; // -4,  3 // 63
    wigner[73] = -temp * pha[1]; // -3,  4 // 73
    wigner[17] = temp * pha[8];  //  4, -3 // 17

    temp = -0.4677071733 * cxp1 * sx3;
    wigner[53] = temp * pha[8];  //  4,  1 // 53
    wigner[27] = -temp * pha[0]; // -4, -1 // 27
    wigner[77] = -temp * pha[5]; //  1,  4 // 77
    wigner[3] = temp * pha[3];   // -1, -4 //  3

    temp = -0.4677071733 * cxm1 * sx3;
    wigner[35] = temp * pha[8];  //  4, -1 // 35
    wigner[45] = -temp * pha[0]; // -4,  1 // 45
    wigner[75] = -temp * pha[3]; // -1,  4 // 75
    wigner[5] = temp * pha[5];   //  1, -4 //  5

    temp = 0.5229125166 * sx3 * sx;
    wigner[44] = temp * pha[8]; //  4,  0 // 44
    wigner[36] = temp * pha[0]; // -4,  0 // 36
    wigner[76] = temp * pha[4]; //  0,  4 // 76
    wigner[4] = temp * pha[4];  //  0, -4 //  4

    temp = -1.4790199458 * sx3 * cx;
    wigner[43] = temp * pha[7];  //  3,  0 // 43
    wigner[37] = -temp * pha[1]; // -3,  0 // 37
    wigner[67] = -temp * pha[4]; //  0,  3 // 67
    wigner[13] = temp * pha[4];  //  0, -3 // 13

    temp = 0.3307189139 * sx2 * cxp12;
    wigner[78] = temp * pha[6]; //  2,  4 // 78
    wigner[2] = temp * pha[2];  // -2, -4 //  2
    wigner[62] = temp * pha[8]; //  4,  2 // 62
    wigner[18] = temp * pha[0]; // -4, -2 // 18

    temp = 0.3307189139 * sx2 * cxm12;
    wigner[6] = temp * pha[6];  //  2, -4 //  6
    wigner[74] = temp * pha[2]; // -2,  4 // 74
    wigner[54] = temp * pha[0]; // -4,  2 // 54
    wigner[26] = temp * pha[8]; //  4, -2 // 26

    temp = 0.4677071733 * cxp12 * sx * (2. * cx - 1.);
    wigner[69] = temp * pha[6];  //  2,  3 // 69
    wigner[11] = -temp * pha[2]; // -2, -3 // 11
    wigner[61] = -temp * pha[7]; //  3,  2 // 61
    wigner[19] = temp * pha[1];  // -3, -2 // 19

    temp = 0.4677071733 * cxm12 * sx * (-2. * cx - 1.);
    wigner[15] = temp * pha[6];  //  2, -3 // 15
    wigner[65] = -temp * pha[2]; // -2,  3 // 65
    wigner[55] = -temp * pha[1]; // -3,  2 // 55
    wigner[25] = temp * pha[7];  //  3, -2 // 25

    temp = 0.25 * cxp12 * (1. - 7. * cxm1 + 7. * cxm12);
    wigner[60] = temp * pha[6]; //  2,  2 // 60
    wigner[20] = temp * pha[2]; // -2, -2 // 20

    temp = 0.25 * cxm12 * (1. - 7. * cxp1 + 7. * cxp12);
    wigner[56] = temp * pha[2]; // -2,  2 // 56
    wigner[24] = temp * pha[6]; //  2, -2 // 24

    temp = 0.3952847075 * sx2 * (7. * cx2 - 1);
    wigner[42] = temp * pha[6]; //  2,  0 // 42
    wigner[38] = temp * pha[2]; // -2,  0 // 38
    wigner[58] = temp * pha[4]; //  0,  2 // 58
    wigner[22] = temp * pha[4]; //  0, -2 // 22

    temp = 0.125 * cxp12 * cxp1 * (-3. + 4. * cx);
    wigner[10] = temp * pha[1]; // -3, -3 // 10
    wigner[70] = temp * pha[7]; //  3,  3 // 70

    temp = 0.125 * cxm12 * cxm1 * (3. + 4. * cx);
    wigner[64] = temp * pha[1]; // -3,  3 // 64
    wigner[16] = temp * pha[7]; //  3, -3 // 16

    temp = 0.3307189139 * cxm1 * cxp12 * (-1. + 4. * cx);
    wigner[12] = temp * pha[3]; // -1, -3 // 12
    wigner[28] = temp * pha[1]; // -3, -1 // 28
    wigner[68] = temp * pha[5]; //  1,  3 // 68
    wigner[52] = temp * pha[7]; //  3,  1 // 52

    temp = 0.3307189139 * cxm12 * cxp1 * (1. + 4. * cx);
    wigner[14] = temp * pha[5]; //  1, -3 // 14
    wigner[46] = temp * pha[1]; // -3,  1 // 46
    wigner[66] = temp * pha[3]; // -1,  3 // 66
    wigner[34] = temp * pha[7]; //  3, -1 // 34

    temp = -0.5590169944 * (4. - 18. * cxm1 + 21. * cxm12 - 7. * cxm12 * cxm1) *
           sx;
    wigner[41] = temp * pha[5];  //  1,  0 // 41
    wigner[39] = -temp * pha[3]; // -1,  0 // 39
    wigner[49] = -temp * pha[4]; //  0,  1 // 49
    wigner[31] = temp * pha[4];  //  0, -1 // 31

    temp = -0.3535533906 * (3. - 10.5 * cxm1 + 7. * cxm12) * sx * cxp1;
    wigner[51] = temp * pha[6];  //  2,  1 // 51
    wigner[29] = -temp * pha[2]; // -2, -1 // 29
    wigner[59] = -temp * pha[5]; //  1,  2 // 59
    wigner[21] = temp * pha[3];  // -1, -2 // 21

    temp = -0.3535533906 * (10. - 17.5 * cxm1 + 7. * cxm12) * sx * cxm1;
    wigner[23] = temp * pha[5];  //  1, -2 // 23
    wigner[57] = -temp * pha[3]; // -1,  2 // 57
    wigner[47] = -temp * pha[2]; // -2,  1 // 47
    wigner[33] = temp * pha[6];  //  2, -1 // 33

    temp = 0.5 * (1. - 9. * cxm1 + 15.75 * cxm12 - 7. * cxm12 * cxm1) * cxp1;
    wigner[30] = temp * pha[3]; // -1, -1 // 30
    wigner[50] = temp * pha[5]; //  1,  1 // 50

    temp = 0.5 * (10. - 30. * cxm1 + 26.25 * cxm12 - 7. * cxm12 * cxm1) * cxm1;
    wigner[32] = temp * pha[5]; //  1, -1 // 32
    wigner[48] = temp * pha[3]; // -1,  1 // 48

    temp = 0.125 * (3. - 30. * cx2 + 35 * cx2 * cx2);
    wigner[40] = temp * pha[4]; //  0,  0 // 40
  }
}

// ❌.. note: (wigner_d) not tested ............................................
double wigner_d(int l, int m1, int m2, double beta) {
  if (l == 2) {
    if (m1 == 2) {
      if (m2 == 2) {
        double cx = cos(beta);
        return ((1 + cx) * (1. + cx) / 4.);
      } else if (m2 == 1) {
        double sx = sin(beta);
        double cx = cos(beta);
        return (-sx * (1. + cx) / 2.);
      } else if (m2 == 0) {
        double sx = sin(beta);
        return (0.6123724355 * sx * sx);
      } else if (m2 == -1) {
        double sx = sin(beta);
        double cx = cos(beta);
        return (-sx * (1. - cx) / 2.);
      } else if (m2 == -2) {
        double cx = cos(beta);
        return ((1 - cx) * (1. - cx) / 4.);
      }
    } else if (m1 == -2) {
      if (m2 == 2) {
        double cx = cos(beta);
        return ((1 - cx) * (1. - cx) / 4.);
      } else if (m2 == 1) {
        double sx = sin(beta);
        double cx = cos(beta);
        return (sx * (1. - cx) / 2.);
      } else if (m2 == 0) {
        double sx = sin(beta);
        return (0.6123724355 * sx * sx);
      } else if (m2 == -1) {
        double sx = sin(beta);
        double cx = cos(beta);
        return (sx * (1. + cx) / 2.);
      } else if (m2 == -2) {
        double cx = cos(beta);
        return ((1 + cx) * (1. + cx) / 4.);
      }
    } else if (m1 == 1) {
      if (m2 == 2) {
        double sx = sin(beta);
        double cx = cos(beta);
        return (sx * (1 + cx) / 2.);
      } else if (m2 == 1) {
        double cx = cos(beta);
        return ((2 * cx * cx + cx - 1.) / 2.);
      } else if (m2 == 0) {
        double sx = sin(beta);
        double cx = cos(beta);
        return (-1.224744871 * sx * cx);
      } else if (m2 == -1) {
        double cx = cos(beta);
        return (-(2 * cx * cx - cx - 1.) / 2.);
      } else if (m2 == -2) {
        double sx = sin(beta);
        double cx = cos(beta);
        return (-sx * (1 - cx) / 2.);
      }
    } else if (m1 == 0) {
      if (m2 == 2) {
        double sx = sin(beta);
        return (0.6123724355 * sx * sx);
      } else if (m2 == 1) {
        double sx = sin(beta);
        double cx = cos(beta);
        return (1.224744871 * sx * cx);
      } else if (m2 == 0) {
        double cx = cos(beta);
        return (1.5 * cx * cx - .5);
      } else if (m2 == -1) {
        double sx = sin(beta);
        double cx = cos(beta);
        return (-1.224744871 * sx * cx);
      } else if (m2 == -2) {
        double sx = sin(beta);
        return (0.6123724355 * sx * sx);
      }
    } else if (m1 == -1) {
      if (m2 == 2) {
        double sx = sin(beta);
        double cx = cos(beta);
        return (sx * (1 - cx) / 2.);
      } else if (m2 == 1) {
        double cx = cos(beta);
        return (-(2 * cx * cx - cx - 1.) / 2.);
      } else if (m2 == 0) {
        double sx = sin(beta);
        double cx = cos(beta);
        return (1.224744871 * sx * cx);
      } else if (m2 == -1) {
        double cx = cos(beta);
        return ((2 * cx * cx + cx - 1.) / 2.);
      } else if (m2 == -2) {
        double sx = sin(beta);
        double cx = cos(beta);
        return (-sx * (1 + cx) / 2.);
      }
    }
  } else {
    double sx = sin(beta / 2.);
    double cx = cos(beta / 2.);
    double sum = 0.;
    int sign = 1;

    for (int k = 0; k <= l - m1; k++) {
      double k1 = (int)(l - m1 - k);
      double k2 = (int)(l + m2 - k);
      double k3 = (int)(k + m1 - m2);

      if (k1 >= 0 && k2 >= 0 && k3 >= 0) {
        int n1 = (int)(2 * l + m2 - m1 - 2 * k);
        int n2 = (int)(m1 - m2 + 2 * k);
        double x = my_power(cx, n1);
        double y = my_power(sx, n2);
        sum += sign * x * y /
               (fac((double)k1) * fac((double)k2) * fac((double)k3) *
                fac((double)k));
      }
      sign = -sign;
    }
    double f = fac(l + m2) * fac(l - m2) * fac(l + m1) * fac(l - m1);
    f = sqrt(f);
    return (sum * f);
  }
  return (0);
}

// ❌.. note: (wigner_d_trig) not tested .......................................
double wigner_d_trig(int l, int m1, int m2, double cx, double sx) {
  if (l == 2) {
    if (m1 == 2) {
      if (m2 == 2) {
        return ((1 + cx) * (1. + cx) / 4.);
      } else if (m2 == 1) {
        return (-sx * (1. + cx) / 2.);
      } else if (m2 == 0) {
        return (0.6123724355 * sx * sx);
      } else if (m2 == -1) {
        return (-sx * (1. - cx) / 2.);
      } else if (m2 == -2) {
        return ((1 - cx) * (1. - cx) / 4.);
      }
    } else if (m1 == -2) {
      if (m2 == 2) {
        return ((1 - cx) * (1. - cx) / 4.);
      } else if (m2 == 1) {
        return (sx * (1. - cx) / 2.);
      } else if (m2 == 0) {
        return (0.6123724355 * sx * sx);
      } else if (m2 == -1) {
        return (sx * (1. + cx) / 2.);
      } else if (m2 == -2) {
        return ((1 + cx) * (1. + cx) / 4.);
      }
    } else if (m1 == 1) {
      if (m2 == 2) {
        return (sx * (1 + cx) / 2.);
      } else if (m2 == 1) {
        return ((2 * cx * cx + cx - 1.) / 2.);
      } else if (m2 == 0) {
        return (-1.224744871 * sx * cx);
      } else if (m2 == -1) {
        return (-(2 * cx * cx - cx - 1.) / 2.);
      } else if (m2 == -2) {
        return (-sx * (1 - cx) / 2.);
      }
    } else if (m1 == 0) {
      if (m2 == 2) {
        return (0.6123724355 * sx * sx);
      } else if (m2 == 1) {
        return (1.224744871 * sx * cx);
      } else if (m2 == 0) {
        return (1.5 * cx * cx - .5);
      } else if (m2 == -1) {
        return (-1.224744871 * sx * cx);
      } else if (m2 == -2) {
        return (0.6123724355 * sx * sx);
      }
    } else if (m1 == -1) {
      if (m2 == 2) {
        return (sx * (1 - cx) / 2.);
      } else if (m2 == 1) {
        return (-(2 * cx * cx - cx - 1.) / 2.);
      } else if (m2 == 0) {
        return (1.224744871 * sx * cx);
      } else if (m2 == -1) {
        return ((2 * cx * cx + cx - 1.) / 2.);
      } else if (m2 == -2) {
        return (-sx * (1 + cx) / 2.);
      }
    }
  }
  return (0);
}

/** ✅ Performs a rank l wigner rotation of the coefficients from the l rank
 * spherical tensors.
 *
 * @param l The rank of the tensor.
 * @param euler_angles A pointer to the array of three euler angles.
 * @param R_in A pointer to the array of coefficients from the l rank tensors of
 *      length 2xl+1 before rotation.
 * @param R_out A pointer to the array of coefficients from the l rank tensors
 *      of length 2xl+1 after rotation.
 */
void single_wigner_rotation(const int l, const double *euler_angles,
                            const complex128 *R_in, complex128 *R_out) {
  int n1 = 2 * l + 1, n2 = n1 * n1, m, mp, k;
  double *wigner = malloc_double(n2);
  // double cos_angle = cos(euler_angles[1]);
  complex128 *temp_initial_vector = malloc_complex128(n1), temp, copy_temp;

  // get wigner matrix corresponding to beta angle
  wigner_d_matrix(l, 1, &euler_angles[1], wigner);

  // copy the initial vector
  cblas_zcopy(n1, R_in, 1, temp_initial_vector, 1);

  // scale the temp initial vector with exp[-I m alpha]
  // orientation at index 0, 1, 2 are alpha, beta, and gamma.
  temp = cos(euler_angles[0]) + I * sin(euler_angles[0]);
  copy_temp = temp;
  for (m = 1; m <= l; m++) {
    temp_initial_vector[l - m] *= temp;
    temp_initial_vector[l + m] *= conj(temp);
    temp *= copy_temp;
  }

  // Apply wigner rotation to the temp inital vector
  k = 0;
  for (m = 0; m < n1; m++) {
    R_out[m] = 0.0;
    for (mp = 0; mp < n1; mp++) {
      R_out[m] += wigner[k++] * temp_initial_vector[mp];
    }
  }

  temp = cos(euler_angles[2]) + I * sin(euler_angles[2]);
  copy_temp = temp;
  for (m = 1; m <= l; m++) {
    R_out[l - m] *= temp;
    R_out[l + m] *= conj(temp);
    temp *= copy_temp;
  }

  free(wigner);
  free(temp_initial_vector);
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
                             const unsigned int n_octants,
                             const double *wigner_2j_matrices,
                             const complex128 *R2,
                             const double *wigner_4j_matrices,
                             const complex128 *R4, complex128 *exp_Im_alpha,
                             complex128 *w2, complex128 *w4) {

  unsigned int j, index_25, index_81, w2_increment, w4_increment = 0;

  complex128 iota = I, negative_iota = -I;
  index_25 = 25 * octant_orientations;
  index_81 = 81 * octant_orientations;
  w2_increment = 5 * octant_orientations;
  if (w4 != NULL) {
    w4_increment = 9 * octant_orientations;
  }

  for (j = 0; j < n_octants; j++) {
    /* Wigner second rank rotation from crystal/common frame to rotor frame.
     */
    __wigner_rotation_2(2, octant_orientations, wigner_2j_matrices,
                        exp_Im_alpha, R2, w2);
    w2 += w2_increment;
    if (sphere) {
      __wigner_rotation_2(2, octant_orientations, &wigner_2j_matrices[index_25],
                          exp_Im_alpha, R2, w2);
      w2 += w2_increment;
    }
    if (w4 != NULL) {
      /* Wigner fourth rank rotation from crystal/common frame to rotor frame.
       */
      __wigner_rotation_2(4, octant_orientations, wigner_4j_matrices,
                          exp_Im_alpha, R4, w4);
      w4 += w4_increment;
      if (sphere) {
        __wigner_rotation_2(4, octant_orientations,
                            &wigner_4j_matrices[index_81], exp_Im_alpha, R4,
                            w4);
        w4 += w4_increment;
      }
    }
    /**
     * Stepping the alpha phase by π/2.
     *
     * The array exp_Im_alpha is a two-dimensional array of shape
     * `4 x number_of_sidebands`, where
     *
     * exp_Im_alpha[m, :] = exp(-I (m-4) alpha) for m=[0, 1, 2, 3]
     * when alpha += π/2
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
    if (!octant) {
      cblas_zscal(octant_orientations, &negative_iota,
                  &exp_Im_alpha[3 * octant_orientations], 1);
      cblas_zdscal(octant_orientations, -1,
                   &exp_Im_alpha[2 * octant_orientations], 1);
      if (w4 != NULL) {
        cblas_zscal(octant_orientations, &iota,
                    &exp_Im_alpha[octant_orientations], 1);
      }
    }
  }
}

/**
 * ✅ Calculates exp(-Im alpha) where alpha is an array of size n.
 * The function accepts cos_alpha = cos(alpha)
 * The result is stored in exp_Im_alpha as m x n matrix where m = [-4,-3,-2,-1]
 */
void get_exp_Im_alpha(const unsigned int n, const double *cos_alpha,
                      const bool allow_fourth_rank, complex128 *exp_Im_alpha) {
  unsigned int i, s_2 = 2 * n, s_3 = 3 * n, s_1 = n;

  // generate sin(alpha) from cos(alpha)
  double *sin_alpha = malloc_double(n);
  for (i = 0; i < n; i++) {
    sin_alpha[i] = sqrt(1.0 - cos_alpha[i] * cos_alpha[i]);
  }

  cblas_dcopy(n, cos_alpha, 1, (double *)&exp_Im_alpha[s_3], 2);
  cblas_dcopy(n, sin_alpha, 1, (double *)&exp_Im_alpha[s_3] + 1, 2);

  vm_double_complex_multiply(n, &exp_Im_alpha[s_3], &exp_Im_alpha[s_3],
                             &exp_Im_alpha[s_2]);

  if (allow_fourth_rank) {
    vm_double_complex_multiply(n, &exp_Im_alpha[s_2], &exp_Im_alpha[s_3],
                               &exp_Im_alpha[s_1]);
    vm_double_complex_multiply(n, &exp_Im_alpha[s_1], &exp_Im_alpha[s_3],
                               exp_Im_alpha);
  }
  free(sin_alpha);
}
