
#include "csa_static_lineshape.h"

static inline void octant(double * spec,
                          double sxx,
                          double syy,
                          double szz,
                          int number_of_points,
                          int nt);

static inline void hemisphere(double * spec,
                          double sxx,
                          double syy,
                          double szz,
                          int number_of_points,
                          int nt);


void lineshape_csa_static(
        double * spec,
        double *cpu_time_,
        int m,
        int nt,
        double fstart,
        double fwidth,
        double iso,
        double aniso,
        double eta,
        int octa,
        int npros) {

  clock_t start, end;
  // Time it
	start = clock();

  double finc, sxx, syy, szz; //, amp1, amp2;
  // int i,j;
  
  // Calculate the principal axis components.
  finc = fwidth/m;
  iso-=fstart;
  iso/=finc;
  aniso/=finc;
  
  sxx = iso -0.5 * aniso* (1.0 + eta);
  syy = iso -0.5 * aniso* (1.0 - eta);
  szz = iso + aniso;

  if (octa){
    octant(spec, sxx, syy, szz, m, nt);
  }
  else{
    hemisphere(spec, sxx, syy, szz, m, nt);
  }
  end = clock();
  cpu_time_[0] = ((double) (end - start)) / (double) CLOCKS_PER_SEC;
}



static inline void octant(double * spec,
                          double sxx,
                          double syy,
                          double szz,
                          int number_of_points,
                          int nt){

  int n_pts = (nt+1)*(nt+2)/2;
  double zero=0.0;


  double* xr = createDouble1DArray(n_pts);
  double* yr = createDouble1DArray(n_pts);
  double* zr = createDouble1DArray(n_pts);
  double* powfreq = createDouble1DArray(n_pts);
  double* amp = createDouble1DArray(n_pts);

  getDirectionCosineSquareOverOctantAndWeights2(nt, xr, yr, zr, amp);

  // calculate the frequencies at various direction cosines. // 
  cblas_daxpy(n_pts, sxx, xr, 1, powfreq, 1);
  cblas_daxpy(n_pts, syy, yr, 1, powfreq, 1);
  cblas_daxpy(n_pts, szz, zr, 1, powfreq, 1);

  // Average over an octant with tenting linear interpolation. //
  powderAverageWithTentingSchemeOverOctant2(
              spec,
              powfreq,
              nt,
              amp,
              &zero,
              number_of_points);

  // free memory
  destroyDouble1DArray(xr);
  destroyDouble1DArray(yr);
  destroyDouble1DArray(zr);
  destroyDouble1DArray(powfreq);
  destroyDouble1DArray(amp);
}


static inline void hemisphere(double * spec,
                        double sxx,
                        double syy,
                        double szz,
                        int number_of_points,
                        int nt){

  int i, j;
  double zero=0.0;
  double** xr = createDouble2DMatrix(nt+1, 2*nt+1);
  double** yr = createDouble2DMatrix(nt+1, 2*nt+1);
  double** zr = createDouble2DMatrix(nt+1, 2*nt+1);
  double** powfreq = createDouble2DMatrix(nt+1, 2*nt+1);
  double** amp = createDouble2DMatrix(nt+1, 2*nt+1);
  getDirectionCosineSquareOverHemishpereAndWeights(nt, xr, yr, zr, amp);
  

  // calculate the frequencies at various direction cosines. // 
  for(i=0; i<nt+1; i++){
    for(j=0; j<nt+1; j++){
      powfreq[i][j] = sxx*xr[i][j] + syy*yr[i][j] + szz*zr[i][j];
    }
  }
  for(i=0; i<nt+1; i++){
    for(j=0; j<2*nt+1; j++){
      powfreq[i][j] = sxx*xr[i][j] + syy*yr[i][j] + szz*zr[i][j];
    }
  }

  // Average over a hemisphere with tenting linear interpolation. //
  powderAverageWithTentingSchemeOverHemisphere(
              spec,
              powfreq,
              nt,
              amp,
              &zero,
              number_of_points);

  // free memory
  destroyDouble2DMatrix(xr);
  destroyDouble2DMatrix(yr);
  destroyDouble2DMatrix(zr);
  destroyDouble2DMatrix(powfreq);
  destroyDouble2DMatrix(amp);
}





