//
//  lebedev131.h
//  Sidebands
//
//  Created by Brennan Walder on 8/8/12.
//  Copyright 2012 Fuck That Shit. All rights reserved.
//

#ifndef Sidebands_lebedev131_h
#define Sidebands_lebedev131_h

void Lebedev_quad_arrays(int N, double alpha[], double beta[], double qweights[]);

void getangle(int N, double *x, double *y, double *z, double *a, double *b);

int Lebedev_Laikov_npoint(int lvalue);

int Lebedev_Laikov_Oh(int n, double a, double b, double v, double *x, double *y, double *z, double *w);

void Lebedev_Laikov_sphere(int N, double *X, double *Y, double *Z, double *W);

#endif
