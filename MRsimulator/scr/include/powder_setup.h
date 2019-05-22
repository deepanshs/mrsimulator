
#include <math.h>
#include <string.h>
#include "c_array.h"
#include <stdio.h>
// #include "core.h"
#include "mkl.h"

extern void getDirectionCosineSquareOverOctantAndWeights2(
          int nt,
          double *xr,
          double *yr,
          double *zr,
          double *amp);

// extern void getDirectionCosineSquareOverOctantAndWeights(
//         int nt,
//         double **xr,
//         double **yr,
//         double **zr,
//         double **rrr);

// extern OCAveragingSchemeDirectionCosineSquareForTenting powderDirectionCosineSquare(int nt);

extern void getDirectionCosineSquareOverHemishpereAndWeights(
        int nt, 
        double ** xr,
        double ** yr,
        double ** zr,
        double ** amp);

// extern void powderDirectionCosineSquare(
//         int nt,
//         double **xr,
//         double **yr,
//         double **zr,
//         double **rrr);

extern void getPolarAngleTrigOverHemisphere(
        int nt,
        double* cosAlpha,
        double* sinAlpha,
        double* cosBeta,
        double* sinBeta,
        double** amp);

void getPolarAngleTrigOverAnOctant(
        int nt,
        double* cosAlpha,
        // double* sinAlpha,
        double* cosBeta,
        // double* sinBeta,
        double* amp);

// void tent(
//         double freq1,
//         double freq2,
//         double freq3,
//         double amp,
//         double *spec,
//         int points,
//         double fstart,
//         double finc);

extern int triangle_interpolation(double *freq1,
          double *freq2,
          double *freq3,
          double *offset,
          double *amp,
          double *spec,
          int *points
          );

// int tent_amp(double *freq1,
//           double *freq2,
//           double *freq3,
//           double *offset,
//           double *amp1,
//           double *amp2,
//           double *amp3,
//           double *spec,
//           int points);

extern void powderAverageWithTentingSchemeOverOctant(
        double *spec,
        double *freq,
        int nt,
        double *amp,
        double *offset,
        int m);

// extern void powderAverageWithTentingSchemeOverOctant(
//         double *spec,
//         double **powfreq,
//         int nt,
//         double **amp,
//         double *offset,
//         int m);

extern void powderAverageWithTentingSchemeOverHemisphere(
        double *spec,
        double **powfreq,
        int nt,
        double **amp,
        double *offset,
        int m);



extern void rasterization(double * grid,
                   double *v0,
                   double *v1,
                   double *v2,
                   int rows,
                   int columns);