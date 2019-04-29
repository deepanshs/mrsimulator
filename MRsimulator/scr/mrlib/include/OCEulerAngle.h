//
//  OCEulerAngle.h
//
//  Created by philip on 4/2/17.
//  Copyright © 2017 Philip Grandinetti. All rights reserved.
//  Contribution: @ 2019 Deepansh J. Srivatava, contact: srivastava.89@osu.edu
//

#ifndef OCEulerAngle_h
#define OCEulerAngle_h

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// OCEuler Type
struct __OCEulerAngle {
    double alphaInRadians;  /* Euler angle alpha */
    double betaInRadians;   /* Euler angle beta */
    double gammaInRadians;  /* Euler angle gamma */
};

typedef struct __OCEulerAngle OCEulerAngle;

void OCEulerAngleShow(OCEulerAngle theAngle);

/*!
 @function OCEulerAngleEqual
 @abstract Returns true if the two angles, 'theAngle1' and 'theAngle2' are equal.
 @param theAngle1 The first angle.
 @param theAngle2 The second angle.
 @result bool.
 */
bool OCEulerAngleEqual(OCEulerAngle theAngle1, OCEulerAngle theAngle2);

/*!
 @function OCEulerAngleZero
 @abstract Returns true if all three angles are zero.
 @param theAngle1 The first angle.
 @result bool.
 */
bool OCEulerAngleZero(OCEulerAngle theAngle);

/*!
 @function OCEulerAngleCreateByAdding
 @abstract Returns the euler angle omega that results from a rotation through the euler theAngle1 and then through the euler theAngle2.
 @param theAngle1 The first angle.
 @param theAngle2 The second angle.
 @result the new angle.
 */
OCEulerAngle OCEulerAngleCreateByAdding(OCEulerAngle theAngle1, OCEulerAngle theAngle2);

OCEulerAngle *OCEulerAngleCreatePowderAngles(char *filename, uint64_t *Nangles, double **weights);

struct __OCHaeberlenConvension {
    double isoInPPM;  /* isotropic chemical shift */
    double zetaInPPM;   /* chemical shift anisotropy */
    double etaInPPM;  /* asymmetry parameter */
};

typedef struct __OCHaeberlenConvension OCHaeberlenConvension;


// OCEuler Type
struct __OCEulerAngleTrig {
    double cos_alphaInRadians;  /* cosine of Euler angle alpha */
    double sin_alphaInRadians;  /* sine of Euler angle alpha */
    double cos_betaInRadians;   /* cosine of Euler angle beta */
    double sin_betaInRadians;   /* sine of Euler angle beta */
    double cos_gammaInRadians;  /* cosine of Euler angle gamma */
    double sin_gammaInRadians;  /* sine of Euler angle gamma */
};
typedef struct __OCEulerAngleTrig OCEulerAngleTrig;

#endif /* OCEulerAngle_h */
