


#ifndef OCPowderScheme_h
#define OCPowderScheme_h

#include "OCEulerAngle.h"

// OCPowderScheme Type
struct __OCPowderScheme {
    // char *name[120];     /* name of the scheme*/
    OCEulerAngle *angles;   /* Euler angles */
    double *weights;        /* Weights of Euler angles */
    int size;               /* number of Euler angles */
};

typedef struct __OCPowderScheme OCPowderScheme;




// OCDirectionCosineSquare Type
struct __OCDirectionCosineSquare {
    double cosSquareX;  /* Euler angle alpha */
    double cosSquareY;   /* Euler angle beta */
    double cosSquareZ;  /* Euler angle gamma */
};
typedef struct __OCDirectionCosineSquare OCDirectionCosineSquare;




// OCPowderSchemeDirectionCosineSquare Type
struct __OCAveragingSchemeDirectionCosineSquareForTenting {
    double **cosSquareX;   /* Euler angles */
    double **cosSquareY;   /* Euler angles */
    double **cosSquareZ;   /* Euler angles */
    double **weights;        /* Weights of Euler angles */
};
typedef struct __OCAveragingSchemeDirectionCosineSquareForTenting \
        OCAveragingSchemeDirectionCosineSquareForTenting;




// OCDirectionCosineSquare Type
struct __OCPolarAngleTrig {
    double cosBeta;  /* cosine of Euler angle beta */
    double sinBeta;   /* sine of Euler angle beta */
    double cosAlpha;  /* cosine of Euler angle alpha */
    double sinAlpha;  /* sine of Euler angle alpha */
};
typedef struct __OCPolarAngleTrig OCPolarAngleTrig;




// OCPowder_direction_cosine_square Type
struct __OCPowderSchemeTrig {
    // char *name[120];     /* name of the scheme*/
    OCPolarAngleTrig **angles;   /* Euler angles */
    double **weights;        /* Weights of Euler angles */
    int size;               /* number of Euler angles */
};
typedef struct __OCPowderSchemeTrig OCPowderSchemeTrig;




OCEulerAngle *OCEulerAngleCreatePowderAngles(char *filename, uint64_t *Nangles, double **weights);

OCPowderScheme OCCreatePowderScheme(int scheme, int size);

#endif /* OCPowderScheme_h */