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
double mypow(double x, int n) {
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

void __wigner_d_matrix(int l, int n, double *angle, double *wigner) {
  vmd_cos(n, angle, angle);
  __wigner_d_matrix_cosine(l, n, angle, wigner);
}

void __wigner_d_matrix_cosine(int l, int n, double *cos_angle, double *wigner) {
  double cx, cx2, sx, temp, *address;
  int i;
  if (l == 2) {
    double t1;
    for (i = 0; i < n; i++) {
      address = &wigner[i * 25];

      cx = cos_angle[i];
      cx2 = cx * cx;
      sx = sqrt(1. - cx2);

      t1 = (1. + cx);
      temp = -sx * t1 / 2.;
      address[19] = temp;  //  2,  1 // 19
      address[5] = -temp;  // -2, -1 //  5
      address[23] = -temp; //  1,  2 // 23
      address[1] = temp;   // -1, -2 //  1

      temp = t1 * t1 / 4.;
      address[24] = temp; //  2,  2 // 24
      address[0] = temp;  // -2, -2 //  0

      t1 = (1. - cx);
      temp = -sx * t1 / 2.;
      address[9] = temp;   //  2, -1 //  9
      address[15] = -temp; // -2,  1 // 15
      address[3] = temp;   //  1, -2 //  3
      address[21] = -temp; // -1,  2 // 21

      temp = t1 * t1 / 4.;
      address[4] = temp;  //  2, -2 //  4
      address[20] = temp; // -2,  2 // 20

      temp = 0.6123724355 * sx * sx;
      address[14] = temp; //  2,  0 // 14
      address[10] = temp; // -2,  0 // 10
      address[22] = temp; //  0,  2 // 22
      address[2] = temp;  //  0, -2 //  2

      temp = 1.224744871 * sx * cx;
      address[13] = -temp; //  1,  0 // 13
      address[17] = temp;  //  0,  1 // 17
      address[7] = -temp;  //  0, -1 //  7
      address[11] = temp;  // -1,  0 // 11

      temp = (2.0 * cx2 + cx - 1.) / 2.;
      address[18] = temp; //  1,  1 // 18
      address[6] = temp;  // -1, -1 //  6

      temp = -(2.0 * cx2 - cx - 1.) / 2.;
      address[8] = temp;  //  1, -1 //  8
      address[16] = temp; // -1,  1 // 16

      address[12] = 1.5 * cx2 - .5; // 0,  0 // 12
    }
  }
  if (l == 4) {
    double sx2, sx3, cxp1, cxm1, cxp12, cxm12, cxm13, cxp13;
    for (i = 0; i < n; i++) {
      address = &wigner[i * 81];

      cx = cos_angle[i];
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
      address[0] = temp;  // -4, -4 //  0
      address[80] = temp; //  4,  4 // 80

      temp = 0.0625 * cxm12 * cxm12;
      address[72] = temp; // -4,  4 // 72
      address[8] = temp;  //  4, -4 //  8

      temp = -0.1767766953 * cxp13 * sx;
      address[1] = temp;   // -3, -4 //  1
      address[9] = -temp;  // -4, -3 //  9
      address[79] = -temp; //  3,  4 // 79
      address[71] = temp;  //  4,  3 //  9

      temp = -0.1767766953 * cxm13 * sx;
      address[7] = temp;   //  3, -4 //  7
      address[63] = -temp; // -4,  3 // 63
      address[73] = -temp; // -3,  4 // 73
      address[17] = temp;  //  4, -3 // 17

      temp = -0.4677071733 * cxp1 * sx3;
      address[53] = temp;  //  4,  1 // 53
      address[27] = -temp; // -4, -1 // 27
      address[77] = -temp; //  1,  4 // 77
      address[3] = temp;   // -1, -4 //  3

      temp = -0.4677071733 * cxm1 * sx3;
      address[35] = temp;  //  4, -1 // 35
      address[45] = -temp; // -4,  1 // 45
      address[75] = -temp; // -1,  4 // 75
      address[5] = temp;   //  1, -4 //  5

      temp = 0.5229125166 * sx3 * sx;
      address[44] = temp; //  4,  0 // 44
      address[36] = temp; // -4,  0 // 36
      address[76] = temp; //  0,  4 // 76
      address[4] = temp;  //  0, -4 //  4

      temp = -1.4790199458 * sx3 * cx;
      address[43] = temp;  //  3,  0 // 43
      address[37] = -temp; // -3,  0 // 37
      address[67] = -temp; //  0,  3 // 67
      address[13] = temp;  //  0, -3 // 13

      temp = 0.3307189139 * sx2 * cxp12;
      address[78] = temp; //  2,  4 // 78
      address[2] = temp;  // -2, -4 //  2
      address[62] = temp; //  4,  2 // 62
      address[18] = temp; // -4, -2 // 18

      temp = 0.3307189139 * sx2 * cxm12;
      address[6] = temp;  //  2, -4 //  6
      address[74] = temp; // -2,  4 // 74
      address[54] = temp; // -4,  2 // 54
      address[26] = temp; //  4, -2 // 26

      temp = 0.4677071733 * cxp12 * sx * (2. * cx - 1.);
      address[69] = temp;  //  2,  3 // 69
      address[11] = -temp; // -2, -3 // 11
      address[61] = -temp; //  3,  2 // 61
      address[19] = temp;  // -3, -2 // 19

      temp = 0.4677071733 * cxm12 * sx * (-2. * cx - 1.);
      address[15] = temp;  //  2, -3 // 15
      address[65] = -temp; // -2,  3 // 65
      address[55] = -temp; // -3,  2 // 55
      address[25] = temp;  //  3, -2 // 25

      temp = 0.25 * cxp12 * (1. - 7. * cxm1 + 7. * cxm12);
      address[60] = temp; //  2,  2 // 60
      address[20] = temp; // -2, -2 // 20

      temp = 0.25 * cxm12 * (1. - 7. * cxp1 + 7. * cxp12);
      address[56] = temp; // -2,  2 // 56
      address[24] = temp; //  2, -2 // 24

      temp = 0.3952847075 * sx2 * (7. * cx2 - 1);
      address[42] = temp; //  2,  0 // 42
      address[38] = temp; // -2,  0 // 38
      address[58] = temp; //  0,  2 // 58
      address[22] = temp; //  0, -2 // 22

      temp = 0.125 * cxp13 * (-3. + 4. * cx);
      address[10] = temp; // -3, -3 // 10
      address[70] = temp; //  3,  3 // 70

      temp = 0.125 * cxm13 * (3. + 4. * cx);
      address[64] = temp; // -3,  3 // 64
      address[16] = temp; //  3, -3 // 16

      temp = 0.3307189139 * cxm1 * cxp12 * (-1. + 4. * cx);
      address[12] = temp; // -1, -3 // 12
      address[28] = temp; // -3, -1 // 28
      address[68] = temp; //  1,  3 // 68
      address[52] = temp; //  3,  1 // 52

      temp = 0.3307189139 * cxm12 * cxp1 * (1. + 4. * cx);
      address[14] = temp; //  1, -3 // 14
      address[46] = temp; // -3,  1 // 46
      address[66] = temp; // -1,  3 // 66
      address[34] = temp; //  3, -1 // 34

      temp = -0.5590169944 * (4. - 18. * cxm1 + 21. * cxm12 - 7. * cxm13) * sx;
      address[41] = temp;  //  1,  0 // 41
      address[39] = -temp; // -1,  0 // 39
      address[49] = -temp; //  0,  1 // 49
      address[31] = temp;  //  0, -1 // 31

      temp = -0.3535533906 * (3. - 10.5 * cxm1 + 7. * cxm12) * sx * cxp1;
      address[51] = temp;  //  2,  1 // 51
      address[29] = -temp; // -2, -1 // 29
      address[59] = -temp; //  1,  2 // 59
      address[21] = temp;  // -1, -2 // 21

      temp = -0.3535533906 * (10. - 17.5 * cxm1 + 7. * cxm12) * sx * cxm1;
      address[23] = temp;  //  1, -2 // 23
      address[57] = -temp; // -1,  2 // 57
      address[47] = -temp; // -2,  1 // 47
      address[33] = temp;  //  2, -1 // 33

      temp = 0.5 * (1. - 9. * cxm1 + 15.75 * cxm12 - 7. * cxm13) * cxp1;
      address[30] = temp; // -1, -1 // 30
      address[50] = temp; //  1,  1 // 50

      temp = 0.5 * (10. - 30. * cxm1 + 26.25 * cxm12 - 7. * cxm13) * cxm1;
      address[32] = temp; //  1, -1 // 32
      address[48] = temp; // -1,  1 // 48

      temp = 0.125 * (3. - 30. * cx2 + 35 * cx2 * cx2);
      address[40] = temp; //  0,  0 // 40
    }
  }
}

void __wigner_rotation(int l, int n, double *wigner, double *cos_alpha,
                       complex128 *R_in, complex128 *R_out) {
  int orientation;
  int n1 = 2 * l + 1, m, i = 0, mp;
  complex128 temp_inital_vector[n1], ph2, pha, ph_temp;
  complex128 *final_vector;
  int n2 = n1 * n1;

  // #pragma omp parallel for
  //   // private(orientation, pha, ph2, m, temp_inital_vector,
  //   //                                  final_vector, i),
  //   //     shared(n, cos_alpha, n1, R_in, l, R_out, wigner, n2)
  for (orientation = 0; orientation < n; orientation++) {

    // calculate the alpha phase, exp(-I m alpha).
    pha.real = cos_alpha[orientation];
    pha.imag = -sqrt(1.0 - cos_alpha[orientation] * cos_alpha[orientation]);

    // note to self: it is important to copy the pha to temporary variable, ph2.
    ph2 = pha;

    // copy the initial vector
    for (m = 0; m < n1; m++) {
      temp_inital_vector[m] = R_in[m];
    }

    // scale the temp initial vector with exp[-I m alpha]
    for (m = 1; m <= l; m++) {
      temp_inital_vector[l + m] = cmult(temp_inital_vector[l + m], ph2);
      ph_temp.real = ph2.real;
      ph_temp.imag = -ph2.imag;
      temp_inital_vector[l - m] = cmult(temp_inital_vector[l - m], ph_temp);
      ph2 = cmult(ph2, pha);
    }

    final_vector = &R_out[orientation * n1];
    i = orientation * n2;
    // Apply wigner rotation to the temp inital vector
    for (m = 0; m < n1; m++) {
      final_vector[m].real = 0.;
      final_vector[m].imag = 0.;
      for (mp = 0; mp < n1; mp++) {
        final_vector[m].real += wigner[i] * temp_inital_vector[mp].real;
        final_vector[m].imag += wigner[i++] * temp_inital_vector[mp].imag;
      }
    }

    // The following blas library is slower than the c for-loops
    // probably because over head in calling blas library is more than the
    // computation of a 5 x 5 matrix-vector multiplication.
    // c for-loops --- 1.19 ms ± 10.4 µs
    // blas        --- 4.19 ms ± 14.8 µs

    // final_vector = &R_out[orientation * n1];
    // cblas_dgemv(CblasRowMajor, CblasNoTrans, n1, n1, 1.0, &wigner[i], n1,
    //             (double *)&temp_inital_vector[0], 2, 0.0,
    //             (double *)&final_vector[0], 2);
    // cblas_dgemv(CblasRowMajor, CblasNoTrans, n1, n1, 1.0, &wigner[i], n1,
    //             (double *)&temp_inital_vector[0] + 1, 2, 0.0,
    //             (double *)&final_vector[0] + 1, 2);
    // i += n2;
  }
}

void __wigner_rotation_2(int l, int n, double *wigner, complex128 *exp_Im_alpha,
                         complex128 *R_in, complex128 *R_out) {
  int orientation;
  int n1 = 2 * l + 1, m, mp, i = 0;
  complex128 *temp_inital_vector = malloc_double_complex(n1);
  complex128 *final_vector, temp, ph_temp;

  // #pragma omp parallel for
  // private(orientation, pha, ph2, m, temp_inital_vector,
  //                                  final_vector, i),
  //     shared(n, cos_alpha, n1, R_in, l, R_out, wigner, n2)
  for (orientation = 0; orientation < n; orientation++) {

    // copy the initial vector
    for (m = 0; m < n1; m++) {
      temp_inital_vector[m] = R_in[m];
    }

    // scale the temp initial vector with exp[-I m alpha]
    for (m = 1; m <= l; m++) {
      temp = exp_Im_alpha[(4 - m) * n + orientation];
      temp_inital_vector[l - m] = cmult(temp_inital_vector[l - m], temp);
      ph_temp.real = temp.real;
      ph_temp.imag = -temp.imag;
      temp_inital_vector[l + m] = cmult(temp_inital_vector[l + m], ph_temp);
    }

    final_vector = &R_out[orientation * n1];
    // Apply wigner rotation to the temp inital vector
    for (m = 0; m < n1; m++) {
      final_vector[m].real = 0.;
      final_vector[m].imag = 0.;
      for (mp = 0; mp < n1; mp++) {
        final_vector[m].real += wigner[i] * temp_inital_vector[mp].real;
        final_vector[m].imag += wigner[i++] * temp_inital_vector[mp].imag;
      }
    }
  }
}

void __wigner_dm0_vector(int l, double beta, double *R_out) {
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

// void full_DLM(complex128 *wigner, int l, double *omega) {
//   if (l == 2) {
//     complex128 pha[5], phg[5];
//     int m;
//     for (m = -2; m <= 2; m++) {
//       pha[m + 2] = cexp(-I * m * omega[0]);
//       // printf("pha %f \n", creal(pha[m+2]));
//       phg[m + 2] = cexp(-I * m * omega[2]);
//     }
//     double cx = cos(omega[1]);
//     double sx = sin(omega[1]);

//     double t1 = (1. + cx);
//     double temp = -sx * t1 / 2.;
//     wigner[19] = temp * pha[4] * phg[3];  //  2,  1  // 19
//     wigner[5] = -temp * pha[0] * phg[1];  // -2, -1 // 5
//     wigner[23] = -temp * pha[3] * phg[4]; //  1,  2 // 23
//     wigner[1] = temp * pha[1] * phg[0];   // -1, -2 // 1

//     temp = t1 * t1 / 4.;
//     wigner[24] = temp * pha[4] * phg[4]; //  2,  2 // 24
//     wigner[0] = temp * pha[0] * phg[0];  // -2, -2 // 0

//     t1 = (1. - cx);
//     temp = -sx * t1 / 2.;
//     wigner[9] = temp * pha[4] * phg[1];   //  2, -1 // 9
//     wigner[15] = -temp * pha[0] * phg[3]; // -2,  1 // 15
//     wigner[3] = temp * pha[3] * phg[0];   //  1, -2 // 3
//     wigner[21] = -temp * pha[1] * phg[4]; // -1,  2 // 21

//     temp = t1 * t1 / 4.;
//     wigner[4] = temp * pha[4] * phg[0];  //  2, -2 // 4
//     wigner[20] = temp * pha[0] * phg[4]; // -2,  2 // 20

//     temp = 0.6123724355 * sx * sx;
//     wigner[14] = temp * pha[4] * phg[2]; //  2,  0 //14
//     wigner[10] = temp * pha[0] * phg[2]; // -2,  0 // 10
//     wigner[22] = temp * pha[2] * phg[4]; //  0,  2 // 22
//     wigner[2] = temp * pha[2] * phg[0];  //  0, -2 // 2

//     temp = 1.224744871 * sx * cx;
//     wigner[13] = -temp * pha[3] * phg[2]; //  1,  0 // 13
//     wigner[17] = temp * pha[2] * phg[3];  //  0,  1 // 17
//     wigner[7] = -temp * pha[2] * phg[1];  //  0, -1 // 7
//     wigner[11] = temp * pha[1] * phg[2];  // -1,  0 // 11

//     temp = (2.0 * cx * cx + cx - 1.) / 2.;
//     wigner[18] = temp * pha[3] * phg[3]; //  1,  1 //18
//     wigner[6] = temp * pha[1] * phg[1];  // -1, -1 // 6

//     temp = -(2.0 * cx * cx - cx - 1.) / 2.;
//     wigner[8] = temp * pha[3] * phg[1];  //  1, -1 // 8
//     wigner[16] = temp * pha[1] * phg[3]; // -1,  1 // 16

//     wigner[12] = 1.5 * cx * cx - .5 * pha[2] * phg[2]; // 0,  0 // 12
//   }
//   if (l == 4) {
//     complex128 pha[9], phg[9];
//     int m;
//     for (m = -4; m <= 4; m++) {
//       pha[m + 4] = cexp(-I * m * omega[0]);
//       // printf("pha %f \n", creal(pha[m+2]));
//       phg[m + 4] = cexp(-I * m * omega[2]);
//     }
//     double cx = cos(omega[1]);
//     double sx = sin(omega[1]);
//     double sx2 = sx * sx, sx3 = sx2 * sx;
//     double cx2 = cx * cx;

//     double cxp1 = (1. + cx), cxm1 = (1. - cx), cxp12 = cxp1 * cxp1,
//            cxm12 = cxm1 * cxm1;

//     // index = (gamma+4)*9 + (alpha+4)
//     double temp = 0.0625 * cxp12 * cxp12;
//     wigner[0] = temp * pha[0];  // -4, -4 //  0
//     wigner[80] = temp * pha[8]; //  4,  4 // 80

//     temp = 0.0625 * cxm12 * cxm12;
//     wigner[72] = temp * pha[0]; // -4,  4 // 72
//     wigner[8] = temp * pha[8];  //  4, -4 //  8

//     temp = -0.1767766953 * cxp12 * cxp1 * sx;
//     wigner[1] = temp * pha[1];   // -3, -4 //  1
//     wigner[9] = -temp * pha[0];  // -4, -3 //  9
//     wigner[79] = -temp * pha[7]; //  3,  4 // 79
//     wigner[71] = temp * pha[8];  //  4,  3 // 71

//     temp = -0.1767766953 * cxm12 * cxm1 * sx;
//     wigner[7] = temp * pha[7];   //  3, -4 //  7
//     wigner[63] = -temp * pha[0]; // -4,  3 // 63
//     wigner[73] = -temp * pha[1]; // -3,  4 // 73
//     wigner[17] = temp * pha[8];  //  4, -3 // 17

//     temp = -0.4677071733 * cxp1 * sx3;
//     wigner[53] = temp * pha[8];  //  4,  1 // 53
//     wigner[27] = -temp * pha[0]; // -4, -1 // 27
//     wigner[77] = -temp * pha[5]; //  1,  4 // 77
//     wigner[3] = temp * pha[3];   // -1, -4 //  3

//     temp = -0.4677071733 * cxm1 * sx3;
//     wigner[35] = temp * pha[8];  //  4, -1 // 35
//     wigner[45] = -temp * pha[0]; // -4,  1 // 45
//     wigner[75] = -temp * pha[3]; // -1,  4 // 75
//     wigner[5] = temp * pha[5];   //  1, -4 //  5

//     temp = 0.5229125166 * sx3 * sx;
//     wigner[44] = temp * pha[8]; //  4,  0 // 44
//     wigner[36] = temp * pha[0]; // -4,  0 // 36
//     wigner[76] = temp * pha[4]; //  0,  4 // 76
//     wigner[4] = temp * pha[4];  //  0, -4 //  4

//     temp = -1.4790199458 * sx3 * cx;
//     wigner[43] = temp * pha[7];  //  3,  0 // 43
//     wigner[37] = -temp * pha[1]; // -3,  0 // 37
//     wigner[67] = -temp * pha[4]; //  0,  3 // 67
//     wigner[13] = temp * pha[4];  //  0, -3 // 13

//     temp = 0.3307189139 * sx2 * cxp12;
//     wigner[78] = temp * pha[6]; //  2,  4 // 78
//     wigner[2] = temp * pha[2];  // -2, -4 //  2
//     wigner[62] = temp * pha[8]; //  4,  2 // 62
//     wigner[18] = temp * pha[0]; // -4, -2 // 18

//     temp = 0.3307189139 * sx2 * cxm12;
//     wigner[6] = temp * pha[6];  //  2, -4 //  6
//     wigner[74] = temp * pha[2]; // -2,  4 // 74
//     wigner[54] = temp * pha[0]; // -4,  2 // 54
//     wigner[26] = temp * pha[8]; //  4, -2 // 26

//     temp = 0.4677071733 * cxp12 * sx * (2. * cx - 1.);
//     wigner[69] = temp * pha[6];  //  2,  3 // 69
//     wigner[11] = -temp * pha[2]; // -2, -3 // 11
//     wigner[61] = -temp * pha[7]; //  3,  2 // 61
//     wigner[19] = temp * pha[1];  // -3, -2 // 19

//     temp = 0.4677071733 * cxm12 * sx * (-2. * cx - 1.);
//     wigner[15] = temp * pha[6];  //  2, -3 // 15
//     wigner[65] = -temp * pha[2]; // -2,  3 // 65
//     wigner[55] = -temp * pha[1]; // -3,  2 // 55
//     wigner[25] = temp * pha[7];  //  3, -2 // 25

//     temp = 0.25 * cxp12 * (1. - 7. * cxm1 + 7. * cxm12);
//     wigner[60] = temp * pha[6]; //  2,  2 // 60
//     wigner[20] = temp * pha[2]; // -2, -2 // 20

//     temp = 0.25 * cxm12 * (1. - 7. * cxp1 + 7. * cxp12);
//     wigner[56] = temp * pha[2]; // -2,  2 // 56
//     wigner[24] = temp * pha[6]; //  2, -2 // 24

//     temp = 0.3952847075 * sx2 * (7. * cx2 - 1);
//     wigner[42] = temp * pha[6]; //  2,  0 // 42
//     wigner[38] = temp * pha[2]; // -2,  0 // 38
//     wigner[58] = temp * pha[4]; //  0,  2 // 58
//     wigner[22] = temp * pha[4]; //  0, -2 // 22

//     temp = 0.125 * cxp12 * cxp1 * (-3. + 4. * cx);
//     wigner[10] = temp * pha[1]; // -3, -3 // 10
//     wigner[70] = temp * pha[7]; //  3,  3 // 70

//     temp = 0.125 * cxm12 * cxm1 * (3. + 4. * cx);
//     wigner[64] = temp * pha[1]; // -3,  3 // 64
//     wigner[16] = temp * pha[7]; //  3, -3 // 16

//     temp = 0.3307189139 * cxm1 * cxp12 * (-1. + 4. * cx);
//     wigner[12] = temp * pha[3]; // -1, -3 // 12
//     wigner[28] = temp * pha[1]; // -3, -1 // 28
//     wigner[68] = temp * pha[5]; //  1,  3 // 68
//     wigner[52] = temp * pha[7]; //  3,  1 // 52

//     temp = 0.3307189139 * cxm12 * cxp1 * (1. + 4. * cx);
//     wigner[14] = temp * pha[5]; //  1, -3 // 14
//     wigner[46] = temp * pha[1]; // -3,  1 // 46
//     wigner[66] = temp * pha[3]; // -1,  3 // 66
//     wigner[34] = temp * pha[7]; //  3, -1 // 34

//     temp = -0.5590169944 * (4. - 18. * cxm1 + 21. * cxm12 - 7. * cxm12 *
//     cxm1) *
//            sx;
//     wigner[41] = temp * pha[5];  //  1,  0 // 41
//     wigner[39] = -temp * pha[3]; // -1,  0 // 39
//     wigner[49] = -temp * pha[4]; //  0,  1 // 49
//     wigner[31] = temp * pha[4];  //  0, -1 // 31

//     temp = -0.3535533906 * (3. - 10.5 * cxm1 + 7. * cxm12) * sx * cxp1;
//     wigner[51] = temp * pha[6];  //  2,  1 // 51
//     wigner[29] = -temp * pha[2]; // -2, -1 // 29
//     wigner[59] = -temp * pha[5]; //  1,  2 // 59
//     wigner[21] = temp * pha[3];  // -1, -2 // 21

//     temp = -0.3535533906 * (10. - 17.5 * cxm1 + 7. * cxm12) * sx * cxm1;
//     wigner[23] = temp * pha[5];  //  1, -2 // 23
//     wigner[57] = -temp * pha[3]; // -1,  2 // 57
//     wigner[47] = -temp * pha[2]; // -2,  1 // 47
//     wigner[33] = temp * pha[6];  //  2, -1 // 33

//     temp = 0.5 * (1. - 9. * cxm1 + 15.75 * cxm12 - 7. * cxm12 * cxm1) * cxp1;
//     wigner[30] = temp * pha[3]; // -1, -1 // 30
//     wigner[50] = temp * pha[5]; //  1,  1 // 50

//     temp = 0.5 * (10. - 30. * cxm1 + 26.25 * cxm12 - 7. * cxm12 * cxm1) *
//     cxm1; wigner[32] = temp * pha[5]; //  1, -1 // 32 wigner[48] = temp *
//     pha[3]; // -1,  1 // 48

//     temp = 0.125 * (3. - 30. * cx2 + 35. * cx2 * cx2);
//     wigner[40] = temp * pha[4]; //  0,  0 // 40
//   }
// }

// // @parameters
// // int l: The angular momentum quantum number
// // complex128 *wigner: a pointer to (2l+1)*(2l+1) size matrix of full
// // Wigner matrix
// void full_DLM_trig(complex128 *wigner, int l, double cosAlpha, double
// sinAlpha,
//                    double cosBeta, double sinBeta) {
//   //    double cosGamma,
//   //    double sinGamma){
//   if (l == 2) {
//     // pha[m+2] holds the value of exp(-I m alpha)
//     // phg[m+2] holds the value of exp(-I m beta)
//     complex128 pha[5]; //, phg[5];
//     pha[2] = 1.0;
//     // phg[2] = 1.0;

//     // for m=-1
//     // ph[1] = cos(-1 alpha) - I sin(-1 alpha)
//     //       = cos(alpha) + I sin(alpha)
//     pha[1] = cosAlpha + I * sinAlpha;
//     // phg[1] = cosGamma + I*sinGamma;

//     // similarly,
//     // for m=1
//     // ph[3] = cos(1 alpha) - I sin(1 alpha)
//     //       = cos(alpha) - I sin(alpha)
//     pha[3] = cosAlpha - I * sinAlpha;
//     // phg[3] = cosGamma - I*sinGamma;

//     // for m= -2
//     // pha[0] = cos(-2 alpha) - I sin(-2 alpha)
//     //        = cos(alpha)^2 - sin(alpha)^2 + I 2 sin(alpha)cos(alpha)
//     //        = ( cos(alpha) + I sin(alpha) )^2
//     //        = (ph[1])^2
//     pha[0] = pha[1] * pha[1];
//     // phg[0] = cpow(phg[1], 2);

//     // similarly for m=2
//     // pha[4] = cos(2 alpha) - I sin(2 alpha)
//     //        = cos(alpha)^2 - sin(alpha)^2 - I 2 sin(alpha)cos(alpha)
//     //        = ( cos(alpha) - I sin(alpha) )^2
//     //        = (ph[3])^2
//     pha[4] = pha[3] * pha[3];
//     // phg[4] = cpow(phg[3], 2);

//     double cx = cosBeta;
//     double sx = sinBeta;
//     double cx2 = cx * cx;
//     // calculating wigner small d terms.

//     double t1 = (1. + cx);
//     double temp = -sx * t1 / 2.;
//     wigner[19] = temp * pha[4];  // * phg[3];    //  2,  1 // 19
//     wigner[5] = -temp * pha[0];  // * phg[1];    // -2, -1 //  5
//     wigner[23] = -temp * pha[3]; // * phg[4];   //  1,  2 // 23
//     wigner[1] = temp * pha[1];   // * phg[0];     // -1, -2 //  1

//     temp = t1 * t1 / 4.;
//     wigner[24] = temp * pha[4]; // * phg[4];    //  2,  2 // 24
//     wigner[0] = temp * pha[0];  // * phg[0];     // -2, -2 //  0

//     t1 = (1. - cx);
//     temp = -sx * t1 / 2.;
//     wigner[9] = temp * pha[4];   // * phg[1];     //  2, -1 //  9
//     wigner[15] = -temp * pha[0]; // * phg[3];   // -2,  1 // 15
//     wigner[3] = temp * pha[3];   // * phg[0];     //  1, -2 //  3
//     wigner[21] = -temp * pha[1]; // * phg[4];   // -1,  2 // 21

//     temp = t1 * t1 / 4.;
//     wigner[4] = temp * pha[4];  // * phg[0];     //  2, -2 //  4
//     wigner[20] = temp * pha[0]; // * phg[4];    // -2,  2 // 20

//     temp = 0.6123724355 * sx * sx;
//     wigner[14] = temp * pha[4]; // * phg[2];    //  2,  0 // 14
//     wigner[10] = temp * pha[0]; // * phg[2];    // -2,  0 // 10
//     wigner[22] = temp * pha[2]; // * phg[4];    //  0,  2 // 22
//     wigner[2] = wigner[22]; // temp * pha[2];// * phg[0];     //  0, -2 //  2

//     temp = 1.224744871 * sx * cx;
//     wigner[13] = -temp * pha[3]; // * phg[2];   //  1,  0 // 13
//     wigner[17] = temp * pha[2];  // * phg[3];    //  0,  1 // 17
//     wigner[7] = -wigner[17];    //-temp * pha[2];// * phg[1];    //  0, -1 //
//     7 wigner[11] = temp * pha[1]; // * phg[2];    // -1,  0 // 11

//     temp = (2.0 * cx2 + cx - 1.) / 2.;
//     wigner[18] = temp * pha[3]; // * phg[3];    //  1,  1 // 18
//     wigner[6] = temp * pha[1];  // * phg[1];     // -1, -1 //  6

//     temp = -(2.0 * cx2 - cx - 1.) / 2.;
//     wigner[8] = temp * pha[3];  // * phg[1];     //  1, -1 //  8
//     wigner[16] = temp * pha[1]; // * phg[3];    // -1,  1 // 16

//     // 0,  0 // 12
//     wigner[12] = 1.5 * cx2 - .5 * pha[2]; // * phg[2];
//   }

//   if (l == 4) {
//     complex128 pha[9]; //, phg[5];
//     pha[4] = 1.0;
//     pha[5] = cosAlpha - I * sinAlpha;
//     pha[3] = cosAlpha + I * sinAlpha;
//     pha[6] = pha[5] * pha[5];
//     pha[2] = conj(pha[6]);
//     pha[7] = pha[6] * pha[5];
//     pha[1] = conj(pha[7]);
//     pha[8] = pha[7] * pha[5];
//     pha[0] = conj(pha[8]);

//     double cx = cosBeta;
//     double sx = sinBeta;
//     double sx2 = sx * sx, sx3 = sx2 * sx;
//     double cx2 = cx * cx;

//     double cxp1 = (1. + cx), cxm1 = (1. - cx), cxp12 = cxp1 * cxp1,
//            cxm12 = cxm1 * cxm1;

//     // index = (gamma+4)*9 + (alpha+4)
//     double temp = 0.0625 * cxp12 * cxp12;
//     wigner[0] = temp * pha[0];  // -4, -4 //  0
//     wigner[80] = temp * pha[8]; //  4,  4 // 80

//     temp = 0.0625 * cxm12 * cxm12;
//     wigner[72] = temp * pha[0]; // -4,  4 // 72
//     wigner[8] = temp * pha[8];  //  4, -4 //  8

//     temp = -0.1767766953 * cxp12 * cxp1 * sx;
//     wigner[1] = temp * pha[1];   // -3, -4 //  1
//     wigner[9] = -temp * pha[0];  // -4, -3 //  9
//     wigner[79] = -temp * pha[7]; //  3,  4 // 79
//     wigner[71] = temp * pha[8];  //  4,  3 // 71

//     temp = -0.1767766953 * cxm12 * cxm1 * sx;
//     wigner[7] = temp * pha[7];   //  3, -4 //  7
//     wigner[63] = -temp * pha[0]; // -4,  3 // 63
//     wigner[73] = -temp * pha[1]; // -3,  4 // 73
//     wigner[17] = temp * pha[8];  //  4, -3 // 17

//     temp = -0.4677071733 * cxp1 * sx3;
//     wigner[53] = temp * pha[8];  //  4,  1 // 53
//     wigner[27] = -temp * pha[0]; // -4, -1 // 27
//     wigner[77] = -temp * pha[5]; //  1,  4 // 77
//     wigner[3] = temp * pha[3];   // -1, -4 //  3

//     temp = -0.4677071733 * cxm1 * sx3;
//     wigner[35] = temp * pha[8];  //  4, -1 // 35
//     wigner[45] = -temp * pha[0]; // -4,  1 // 45
//     wigner[75] = -temp * pha[3]; // -1,  4 // 75
//     wigner[5] = temp * pha[5];   //  1, -4 //  5

//     temp = 0.5229125166 * sx3 * sx;
//     wigner[44] = temp * pha[8]; //  4,  0 // 44
//     wigner[36] = temp * pha[0]; // -4,  0 // 36
//     wigner[76] = temp * pha[4]; //  0,  4 // 76
//     wigner[4] = temp * pha[4];  //  0, -4 //  4

//     temp = -1.4790199458 * sx3 * cx;
//     wigner[43] = temp * pha[7];  //  3,  0 // 43
//     wigner[37] = -temp * pha[1]; // -3,  0 // 37
//     wigner[67] = -temp * pha[4]; //  0,  3 // 67
//     wigner[13] = temp * pha[4];  //  0, -3 // 13

//     temp = 0.3307189139 * sx2 * cxp12;
//     wigner[78] = temp * pha[6]; //  2,  4 // 78
//     wigner[2] = temp * pha[2];  // -2, -4 //  2
//     wigner[62] = temp * pha[8]; //  4,  2 // 62
//     wigner[18] = temp * pha[0]; // -4, -2 // 18

//     temp = 0.3307189139 * sx2 * cxm12;
//     wigner[6] = temp * pha[6];  //  2, -4 //  6
//     wigner[74] = temp * pha[2]; // -2,  4 // 74
//     wigner[54] = temp * pha[0]; // -4,  2 // 54
//     wigner[26] = temp * pha[8]; //  4, -2 // 26

//     temp = 0.4677071733 * cxp12 * sx * (2. * cx - 1.);
//     wigner[69] = temp * pha[6];  //  2,  3 // 69
//     wigner[11] = -temp * pha[2]; // -2, -3 // 11
//     wigner[61] = -temp * pha[7]; //  3,  2 // 61
//     wigner[19] = temp * pha[1];  // -3, -2 // 19

//     temp = 0.4677071733 * cxm12 * sx * (-2. * cx - 1.);
//     wigner[15] = temp * pha[6];  //  2, -3 // 15
//     wigner[65] = -temp * pha[2]; // -2,  3 // 65
//     wigner[55] = -temp * pha[1]; // -3,  2 // 55
//     wigner[25] = temp * pha[7];  //  3, -2 // 25

//     temp = 0.25 * cxp12 * (1. - 7. * cxm1 + 7. * cxm12);
//     wigner[60] = temp * pha[6]; //  2,  2 // 60
//     wigner[20] = temp * pha[2]; // -2, -2 // 20

//     temp = 0.25 * cxm12 * (1. - 7. * cxp1 + 7. * cxp12);
//     wigner[56] = temp * pha[2]; // -2,  2 // 56
//     wigner[24] = temp * pha[6]; //  2, -2 // 24

//     temp = 0.3952847075 * sx2 * (7. * cx2 - 1);
//     wigner[42] = temp * pha[6]; //  2,  0 // 42
//     wigner[38] = temp * pha[2]; // -2,  0 // 38
//     wigner[58] = temp * pha[4]; //  0,  2 // 58
//     wigner[22] = temp * pha[4]; //  0, -2 // 22

//     temp = 0.125 * cxp12 * cxp1 * (-3. + 4. * cx);
//     wigner[10] = temp * pha[1]; // -3, -3 // 10
//     wigner[70] = temp * pha[7]; //  3,  3 // 70

//     temp = 0.125 * cxm12 * cxm1 * (3. + 4. * cx);
//     wigner[64] = temp * pha[1]; // -3,  3 // 64
//     wigner[16] = temp * pha[7]; //  3, -3 // 16

//     temp = 0.3307189139 * cxm1 * cxp12 * (-1. + 4. * cx);
//     wigner[12] = temp * pha[3]; // -1, -3 // 12
//     wigner[28] = temp * pha[1]; // -3, -1 // 28
//     wigner[68] = temp * pha[5]; //  1,  3 // 68
//     wigner[52] = temp * pha[7]; //  3,  1 // 52

//     temp = 0.3307189139 * cxm12 * cxp1 * (1. + 4. * cx);
//     wigner[14] = temp * pha[5]; //  1, -3 // 14
//     wigner[46] = temp * pha[1]; // -3,  1 // 46
//     wigner[66] = temp * pha[3]; // -1,  3 // 66
//     wigner[34] = temp * pha[7]; //  3, -1 // 34

//     temp = -0.5590169944 * (4. - 18. * cxm1 + 21. * cxm12 - 7. * cxm12 *
//     cxm1) *
//            sx;
//     wigner[41] = temp * pha[5];  //  1,  0 // 41
//     wigner[39] = -temp * pha[3]; // -1,  0 // 39
//     wigner[49] = -temp * pha[4]; //  0,  1 // 49
//     wigner[31] = temp * pha[4];  //  0, -1 // 31

//     temp = -0.3535533906 * (3. - 10.5 * cxm1 + 7. * cxm12) * sx * cxp1;
//     wigner[51] = temp * pha[6];  //  2,  1 // 51
//     wigner[29] = -temp * pha[2]; // -2, -1 // 29
//     wigner[59] = -temp * pha[5]; //  1,  2 // 59
//     wigner[21] = temp * pha[3];  // -1, -2 // 21

//     temp = -0.3535533906 * (10. - 17.5 * cxm1 + 7. * cxm12) * sx * cxm1;
//     wigner[23] = temp * pha[5];  //  1, -2 // 23
//     wigner[57] = -temp * pha[3]; // -1,  2 // 57
//     wigner[47] = -temp * pha[2]; // -2,  1 // 47
//     wigner[33] = temp * pha[6];  //  2, -1 // 33

//     temp = 0.5 * (1. - 9. * cxm1 + 15.75 * cxm12 - 7. * cxm12 * cxm1) * cxp1;
//     wigner[30] = temp * pha[3]; // -1, -1 // 30
//     wigner[50] = temp * pha[5]; //  1,  1 // 50

//     temp = 0.5 * (10. - 30. * cxm1 + 26.25 * cxm12 - 7. * cxm12 * cxm1) *
//     cxm1; wigner[32] = temp * pha[5]; //  1, -1 // 32 wigner[48] = temp *
//     pha[3]; // -1,  1 // 48

//     temp = 0.125 * (3. - 30. * cx2 + 35 * cx2 * cx2);
//     wigner[40] = temp * pha[4]; //  0,  0 // 40
//   }
// }

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
        double x = mypow(cx, n1);
        double y = mypow(sx, n2);
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

// cx = cos(beta) and sx = sin(beta)
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
