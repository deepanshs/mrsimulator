#include "c_array.h" 
#include <complex.h>
#include "OCPowderScheme.h"

// float 1d array
float* createFloat1DArray(int m)
{
    float* values = calloc(m, sizeof(float));
    return values;
}

void destroyFloat1DArray(float* arr)
{
    free(arr);
}


// complex float 1d array
float complex* createFloatComplex1DArray(int m)
{
    float complex* values = calloc(m, sizeof(float complex));
    return values;
}

void destroyFloatComplex1DArray(float complex* arr)
{
    free(arr);
}


// double 1d array
double* createDouble1DArray(int m)
{
    double* values = calloc(m, sizeof(double));
    return values;
}

void destroyDouble1DArray(double* arr)
{
    free(arr);
}


// complex double 1d array
double complex* createDoubleComplex1DArray(int m)
{
    double complex* values = calloc(m, sizeof(double complex));
    return values;
}

void destroyDoubleComplex1DArray(double complex* arr)
{
    free(arr);
}



// float 2d matrix
float** createFloat2DMatrix(int n, int m)
{
    float* values = calloc(m*n, sizeof(float));
    float** rows = malloc(n*sizeof(float*));
    for (int i=0; i<n; ++i)
    {
        rows[i] = values + i*m;
    }
    return rows;
}

void destroyFloat2DMatrix(float** arr)
{
    free(*arr);
    free(arr);
}


// double array
double** createDouble2DMatrix(int n, int m)
{
  double* values = calloc(m*n, sizeof(double));
    double** rows = malloc(n*sizeof(double*));
    for (int i=0; i<n; ++i)
    {
        rows[i] = values + i*m;
    }
    return rows;
}

void destroyDouble2DMatrix(double** arr)
{
    free(*arr);
    free(arr);
}


// double 3d array
double*** createDouble3DArray(int n, int m, int o)
{
//   double* values = calloc(n*m*o, sizeof(double));
//     double ** columns = malloc(m*sizeof(double*));
//         double*** rows = malloc(n*sizeof(double**));
//         for (int i=0; i<m; ++i)
//         {
//             columns[i] = values + i*m;
//             for(int j=0;  j<n;  j++){
//                 rows[i][j] = columns[i] + j*n
//             }
//         }
//         return rows;
// }
    double *** p= malloc(n * sizeof(double ***));
    for(int i=0;  i<n;  i++)
        {
            p[i]=(double **)malloc(m * sizeof(double *));
            
            for(int j=0;  j<m;  j++)
            {
                p[i][j]=(double *)malloc(o * sizeof(double));
            }
        }
    return p;
}

void destroyDouble3DArray(double*** arr)
{
    free(**arr);
    free(*arr);
    free(arr);
}

// Create 2D array of polar angle trigs
OCPolarAngleTrig** create2DOCPolarAngleTrigArray(int n, int m)
{
    OCPolarAngleTrig* values = calloc(m*n, sizeof(OCPolarAngleTrig));
    OCPolarAngleTrig** rows = malloc(n*sizeof(OCPolarAngleTrig*));
    for (int i=0; i<n; ++i)
    {
        rows[i] = values + i*m;
    }
    return rows;
}

void destroy2DOCPolarAngleTrigArray(OCPolarAngleTrig** arr)
{
    free(*arr);
    free(arr);
}


// Create 2D array of directon cosines
OCDirectionCosineSquare** create2DOCDirectionCosineSquareArray(int n, int m)
{
    OCDirectionCosineSquare* values = calloc(m*n, sizeof(OCDirectionCosineSquare));
    OCDirectionCosineSquare** rows = malloc(n*sizeof(OCDirectionCosineSquare*));
    for (int i=0; i<n; ++i)
    {
        rows[i] = values + i*m;
    }
    return rows;
}

void destroy2DOCDirectionCosineSquareArray(OCDirectionCosineSquare** arr)
{
    free(*arr);
    free(arr);
}

// // main 
// int main() {
//   int i, j, m, n;
//   m=13;
//   n=2;
//   double** arr = createDouble2DMatrix(n,m);
//   // arr = array(3,3);
//   for ( i = 0; i < n; i++ ) {
//     for ( j = 0; j < m; j++ ) {
//       arr[i][j] = 1.0;
//       printf("%f\n", arr[i][j]);
//     }
//   }
//   printf("\n new \n");
//   modify2DMatrix(arr, n, m);

//   for ( i = 0; i < n; i++ ) {
//     for ( j = 0; j < m; j++ ) {
//       printf("%f\n", arr[i][j]);
//     }
//   }
//   destroyDouble2DMatrix(arr);



//   double* arr1 = createDouble1DArray(m);
//   // arr = array(3,3);
//   for ( i = 0; i < m; i++ ) {
//       arr1[i] = 1.0;
//       printf("%f\n", arr1[i]);
//     }
  
//   printf("\n new \n");
//   modify1DArray(arr1, m);

//   for ( i = 0; i < m; i++ ) {
//       printf("%f\n", arr1[i]);
//     }

//   destroyDouble1DArray(arr1);

//   return 0;
// }


// test modify array
void modify2DMatrix(double **arr, int n, int m)
{
  int i, j;
  for( i = 0; i < n; i++) {
    for( j = 0; j < m; j++) {
      arr[i][j] = i+j;
    }
  }
}

void modify1DArray(double *arr, int m)
{
  int i;
  for( i = 0; i < m; i++) {
      arr[i] = i;
    }
}