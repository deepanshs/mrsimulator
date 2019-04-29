//
//  OCEulerAngle.c
//
//  Created by philip on 4/2/17.
//  Copyright © 2017 Philip Grandinetti. All rights reserved.
//  Contribution: Deepansh J. Srivatava. contact: srivastava.89@osu.edu
//

#include "OCEulerAngle.h"

bool OCEulerAngleEqual(OCEulerAngle theAngle1, OCEulerAngle theAngle2)
{
    if(theAngle1.alphaInRadians != theAngle2.alphaInRadians) return false;
    if(theAngle1.betaInRadians != theAngle2.betaInRadians) return false;
    if(theAngle1.gammaInRadians != theAngle2.gammaInRadians) return false;
    return true;
}

// OCEulerAngle OCEulerAngleCreate(double alphaInRadians, double betaInRadians, double gammaInRadians)
// {
//     return (OCEulerAngle) {alphaInRadians,betaInRadians,gammaInRadians};
// }

bool OCEulerAngleZero(OCEulerAngle theAngle)
{
    bool value = false;
    if(theAngle.alphaInRadians == 0.0 && theAngle.betaInRadians == 0.0 &&  theAngle.gammaInRadians == 0.0) value=true;
    return value;
}

void OCEulerAngleShow(OCEulerAngle theAngle)
{
    printf("Euler angle with alpha=%lf, beta=%lf, gamma=%lf\n",theAngle.alphaInRadians, theAngle.betaInRadians,theAngle.gammaInRadians);
}

OCEulerAngle OCEulerAngleCreateByAdding(OCEulerAngle theAngle1, OCEulerAngle theAngle2)
{
    /*	Calculate some useful numbers */
    double cosb1 = cos(theAngle1.betaInRadians);
    double cosb2 = cos(theAngle2.betaInRadians);
    double sinb1 = sin(theAngle1.betaInRadians);
    double sinb2 = sin(theAngle2.betaInRadians);
    double sina1g2 = sin(theAngle1.alphaInRadians + theAngle2.gammaInRadians);
    double cosa1g2 = cos(theAngle1.alphaInRadians + theAngle2.gammaInRadians);
    
    double alpha = 0;
    double beta = acos(cosb1*cosb2-sinb1*sinb2*cosa1g2);
    double gamma = 0;
    
    double sinb = sin(beta);
    
    if(sinb != 0.) {
        double cosaa2 = cosa1g2 * sinb1*cosb2/sinb + sinb2*cosb1/sinb;
        double sinaa2 = sina1g2 * sinb1/sinb;
        double aa2 = aa2 = acos(cosaa2);
        if(cosaa2>1.0) aa2 = 0.;
        else if(cosaa2<-1.0) aa2 = M_PI;
        if(sinaa2<0.) aa2 = 2.*M_PI - aa2;
        alpha = theAngle2.alphaInRadians + aa2;
        
        double cosgg1 = cosa1g2 * sinb2*cosb1/sinb + sinb1*cosb2/sinb;
        double singg1 = sina1g2 * sinb2/sinb;
        double gg1 = acos(cosgg1);
        if(cosgg1>1.0) gg1 = 0.;
        else if(cosgg1<-1.0) gg1 = M_PI;
        if(singg1<0.) gg1 = 2.*M_PI - gg1;
        gamma = theAngle1.gammaInRadians + gg1;
    }
    else {
        /* This isn't right */
        alpha = theAngle2.alphaInRadians;
        gamma = theAngle1.gammaInRadians;
    }
    
    OCEulerAngle result = {alpha,beta,gamma};
    return result;
}


