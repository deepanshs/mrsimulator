
#include "powder_setup.h"
#include "OCPowderScheme.h"

void getDirectionCosineSquareOverOctantAndWeights2(
          int nt,
          double *xr,
          double *yr,
          double *zr,
          double *amp) {

  int i, j, k=0;
  double x2, y2, z2, r2, scale = 1.0; // /((double)(nt+1)*(nt+2)/2);

  /* Do the (x + y + z = nt) face of the octahedron
  !z -> 0 to nt-1
  !y -> 0 to nt-z
  !x -> nt - y - z
  !*/

  for( j = 0; j <= nt-1; j++) {
    for( i = 0; i <= nt-j; i++) {
      // x = nt-i-j;
      // y = i;
      // z = j;
      x2 = pow(nt-i-j, 2);
      y2 = pow(i, 2);
      z2 = pow(j, 2);
      r2 = x2 + y2 + z2;
      xr[k] = x2/r2;
      yr[k] = y2/r2;
      zr[k] = z2/r2;
      amp[k] = scale/(r2*sqrt(r2));
      k++;
    }
  }

  xr[k] = 0.0;
  yr[k] = 0.0;
  zr[k] = 1.0;
  r2 = nt;
  amp[k] = scale/(r2*r2*r2);
}


// void getDirectionCosineSquareOverOctantAndWeights(
//           int nt,
//           double **xr,
//           double **yr,
//           double **zr,
//           double **amp) {

//   int i, j;
//   double x2, y2, z2, r2;

//   /* Do the (x + y + z = nt) face of the octahedron
//   !z -> 0 to nt-1
//   !y -> 0 to nt-z
//   !x -> nt - y - z
//   !*/

//   for( j = 0; j <= nt-1; j++) {
//     for( i = 0; i <= nt-j; i++) {
//       // x = nt-i-j;
//       // y = i;
//       // z = j;
//       x2 = pow(nt-i-j, 2);
//       y2 = pow(i, 2);
//       z2 = pow(j, 2);
//       r2 = x2 + y2 + z2;
//       xr[i][j] = x2/r2;
//       yr[i][j] = y2/r2;
//       zr[i][j] = z2/r2;
//       amp[i][j] = 1.0/(r2*sqrt(r2));
//     }
//   }

//   xr[0][nt] = 0.0;
//   yr[0][nt] = 0.0;
//   zr[0][nt] = 1.0;
//   r2 = nt;
//   amp[0][nt] = 1.0/(r2*r2*r2);
// }

void getDirectionCosineSquareOverHemishpereAndWeights(
            int nt, 
            double ** xr,
            double ** yr,
            double ** zr,
            double ** amp){

  int i, j;
  double x2, y2, z2, r2;

  /* Do the (x + y + z = nt) face of the octahedron
  !z -> 0 to nt-1
  !y -> 0 to nt-z
  !x -> nt - y - z
  !*/
  
  for( j = 0; j <= nt-1; j++) {
    for( i = 0; i <= nt-j; i++) {
      // x = nt-i-j;
      // y = i;
      // z = j;
      x2 = pow(nt-i-j, 2);
      y2 = pow(i, 2);
      z2 = pow(j, 2);
      r2 = x2 + y2 + z2;
      xr[i][j] = x2/r2;
      yr[i][j] = y2/r2;
      zr[i][j] = z2/r2;
      amp[i][j] = 1.0/(r2*sqrt(r2));
    }

  /* Do the (-x + y + z = nt) face of the octahedron
  !z -> 0 to nt-1
  !y -> 0 to nt-z
  !x -> -nt + y + z
  !*/
    for( i = nt-j+1; i <= nt; i++) {
      // x = nt-i-j;
      // y = nt-j;
      // z = nt-i;
      x2 = pow(nt-i-j,2);
      y2 = pow(nt-j, 2);
      z2 = pow(nt-i, 2);
      r2 = x2 + y2 + z2;
      xr[i][j] = x2/r2;
      yr[i][j] = y2/r2;
      zr[i][j] = z2/r2;
      amp[i][j] = 1.0/(r2*sqrt(r2));
    }
  }

  /* Do the (-x - y + z = nt) face of the octahedron
  !*/


  for( j = nt; j < 2*nt; j++) {
    for( i = j-nt+1; i < nt; i++) {
      // x = -nt-i+j;
      // y = nt-j;
      // z = nt-i;
      x2 = pow(-nt-i+j, 2); // x*x;
      y2 = pow(nt-j, 2); // y*y;
      z2 = pow(nt-i, 2); //z*z;
      r2 = x2 + y2 + z2;
      xr[i][j] = x2/r2;
      yr[i][j] = y2/r2;
      zr[i][j] = z2/r2;
      amp[i][j] = 1.0/(r2*sqrt(r2));
    }

  /* Do the (x - y + z = nt) face of the octahedron
  !*/

    for( i = 1; i <= j-nt; i++) {
      // x = -nt-i+j;
      // y = -i;
      // z = 2*nt-j;
      x2 = pow(-nt-i+j, 2); //x*x;
      y2 = pow(-i, 2); //y*y;
      z2 = pow(2*nt-j, 2); //z*z;
      r2 = x2 + y2 + z2;
      xr[i][j] = x2/r2;
      yr[i][j] = y2/r2;
      zr[i][j] = z2/r2;
      amp[i][j] = 1.0/(r2*sqrt(r2));
    }
  }

  xr[0][nt] = 0.0;
  yr[0][nt] = 0.0;
  zr[0][nt] = 1.0;
  r2 = nt;
  amp[0][nt] = 1.0/(r2*r2*r2);

  for( j = 0; j < nt; j++){
      i = 2*nt-j;
      xr[0][i] = xr[0][j];
      yr[0][i] = yr[0][j];
      zr[0][i] = zr[0][j];
      amp[0][i] = amp[0][j];
  }

  for( i = 0; i <= nt; i++) {
      xr[nt][nt+i] = xr[i][0];
      yr[nt][nt+i] = yr[i][0];
      zr[nt][nt+i] = zr[i][0];
      amp[nt][nt+i] = amp[i][0];
  }

  i = 2*nt;
  for( j = 1; j < nt; j++) {
      xr[nt-j][i] = xr[nt][j];
      yr[nt-j][i] = yr[nt][j];
      zr[nt-j][i] = zr[nt][j];
      amp[nt-j][i] = amp[nt][j];
  }


//   if (accurateArea == 1){

//     double x12, y12, z12, x13, y13, z13, pyz, pzx, pxy, temp0, temp;
//     int k=0;

//     // double **x = createDouble2DMatrix( (nt+1), (2*nt+1) );
//     // double **y = createDouble2DMatrix( (nt+1), (2*nt+1) );
//     // double **z = createDouble2DMatrix( (nt+1), (2*nt+1) );

//     // vdSqrt(points, &zr[0][0], &z[0][0]);

//     // // Evaluate sqrt of xr
//     // vdSqrt(points, &xr[0][0], &x[0][0]);

//     // // Evaluate sqrt of xr
//     // vdSqrt(points, &yr[0][0], &y[0][0]);

//     temp0 = amp[0][0] + amp[0][1] + amp[1][0];
//     for (i=0; i<=nt-1; i++){
//       for (j=0; j<=nt-1; j++){

//         // First
//         x12 = sqrt(xr[i][j+1]) - sqrt(xr[i+1][j]);
//         y12 = sqrt(yr[i][j+1]) - sqrt(yr[i+1][j]);
//         z12 = sqrt(zr[i][j+1]) - sqrt(zr[i+1][j]);

//         x13 = sqrt(xr[i][j+1]) - sqrt(xr[i][j]);
//         y13 = sqrt(yr[i][j+1]) - sqrt(yr[i][j]);
//         z13 = sqrt(zr[i][j+1]) - sqrt(zr[i][j]);

//         pyz = y12*z13-z12*y13;
//         pzx = z12*x13-x12*z13;
//         pxy = x12*y13-y12*x13;

//         amp1[k] = sqrt(pyz*pyz + pzx*pzx + pxy*pxy); // * 0.5
//         temp = amp[i][j] + amp[i][j+1] + amp[i+1][j];
//         printf("%f %f %f \n", amp1[k]/amp1[0], temp/temp0, (amp1[k]/amp1[0]) - (temp/temp0) );
//         k++;


//         // Second
//         x13 = sqrt(xr[i][j+1]) - sqrt(xr[i+1][j+1]);
//         y13 = sqrt(yr[i][j+1]) - sqrt(yr[i+1][j+1]);
//         z13 = sqrt(zr[i][j+1]) - sqrt(zr[i+1][j+1]);

//         // x13 = xr[i][j+1] - xr[i+1][j];
//         // y13 = yr[i][j+1] - yr[i+1][j];
//         // z13 = zr[i][j+1] - zr[i+1][j];

//         pyz = y12*z13-z12*y13;
//         pzx = z12*x13-x12*z13;
//         pxy = x12*y13-y12*x13;

//         amp1[k] = sqrt(pyz*pyz + pzx*pzx + pxy*pxy); // * 0.5
//         temp = amp[i][j+1] + amp[i+1][j+1] + amp[i+1][j];
//         printf("%f %f %f \n", amp1[k]/amp1[0], temp/temp0, (amp1[k]/amp1[0]) - (temp/temp0) );
//         k++;
//       }
//     }

//     for (i=0; i<=nt-1; i++){
//       for (j=nt; j<=2*nt-1; j++){

//         // Third
//         x12 = sqrt(xr[i][j]) - sqrt(xr[i+1][j+1]);
//         y12 = sqrt(yr[i][j]) - sqrt(yr[i+1][j+1]);
//         z12 = sqrt(zr[i][j]) - sqrt(zr[i+1][j+1]);

//         x13 = sqrt(xr[i][j]) - sqrt(xr[i+1][j]);
//         y13 = sqrt(yr[i][j]) - sqrt(yr[i+1][j]);
//         z13 = sqrt(zr[i][j]) - sqrt(zr[i+1][j]);

//         pyz = y12*z13-z12*y13;
//         pzx = z12*x13-x12*z13;
//         pxy = x12*y13-y12*x13;

//         amp1[k] = sqrt(pyz*pyz + pzx*pzx + pxy*pxy); // * 0.5
//         temp = amp[i][j] + amp[i+1][j+1] + amp[i+1][j];
//         printf("%f %f %f \n", amp1[k]/amp1[0], temp/temp0, (amp1[k]/amp1[0]) - (temp/temp0) );
//         k++;

//         // Fourth
//         // x12 = xr[i][j] - xr[i+1][j+1];
//         // y12 = yr[i][j] - yr[i+1][j+1];
//         // z12 = zr[i][j] - zr[i+1][j+1];

//         x13 = sqrt(xr[i][j]) - sqrt(xr[i][j+1]);
//         y13 = sqrt(yr[i][j]) - sqrt(yr[i][j+1]);
//         z13 = sqrt(zr[i][j]) - sqrt(zr[i][j+1]);

//         pyz = y12*z13-z12*y13;
//         pzx = z12*x13-x12*z13;
//         pxy = x12*y13-y12*x13;

//         amp1[k] = sqrt(pyz*pyz + pzx*pzx + pxy*pxy); // * 0.5
//         temp = amp[i][j] + amp[i+1][j+1] + amp[i][j+1];
//         printf("%f %f %f \n", amp1[k]/amp1[0], temp/temp0, (amp1[k]/amp1[0]) - (temp/temp0) );
//         k++;
//       }
//     }
//   }
//   else{
//     int k=0;
//     for (i=0; i<=nt-1; i++){
//       for (j=0; j<=nt-1; j++){

//         // First
//         amp1[k++] = amp[i][j] + amp[i][j+1] + amp[i+1][j]; // * 0.5

//         // Second
//         amp1[k++]= amp[i][j+1] + amp[i+1][j] + amp[i+1][j+1]; // * 0.5
//       }
//     }
//     for (i=0; i<=nt-1; i++){
//       for (j=nt; j<=2*nt-1; j++){

//         // Third
//         amp1[k++] = amp[i][j] + amp[i+1][j+1] + amp[i+1][j]; // * 0.5

//         // Fourth
//         amp1[k++] = amp[i][j] + amp[i+1][j+1] + amp[i][j+1]; // * 0.5
//       }
//     }
//   }
//   destroyDouble2DMatrix(amp);
}

void getPolarAngleTrigOverAnOctant(
        int nt,
        double* cosAlpha,
        double* sinAlpha,
        double* cosBeta,
        double* sinBeta,
        double* amp)
{

  // int i, j;

  int points = (nt+1) * (nt+2)/2;
  double* xr = createDouble1DArray( points );
  double* yr = createDouble1DArray( points );
  double* zr = createDouble1DArray( points );
  
  getDirectionCosineSquareOverOctantAndWeights2(nt, xr, yr, zr, amp);

  // // OCPolarAngleTrig **trigs = create2DOCPolarAngleTrigArray(nt+1, 2*nt+1);

  // Evaluate sqrt of zr to get cos(beta)
  vdSqrt(points, &zr[0], &cosBeta[0]);

   // exaluate A = x + y
  vdAdd(points, &xr[0], &yr[0], &sinBeta[0]);
  // Take sqrt of A to get sin(beta)
  vdSqrt(points, &sinBeta[0], &sinBeta[0]);

  // Evaluate sqrt of xr
  vdSqrt(points, &xr[0], &xr[0]);

  // Evaluate sqrt of xr
  vdSqrt(points, &yr[0], &yr[0]);

  vdDiv(points-1, xr, sinBeta, cosAlpha );
  vdDiv(points-1, yr, sinBeta, sinAlpha );

  cosAlpha[points-1] = 1.0;
  sinAlpha[points-1] = 0.0;

  // int ii=0;
  // for( i = 0; i < points; i++) {
  //   if (sinBeta[i] != 0.){
  //     cosAlpha[i] = xr[i] / sinBeta[i];
  //     sinAlpha[i] = yr[i] / sinBeta[i];
  //   }
  //   else{
  //     cosAlpha[i] = 1.0;
  //     sinAlpha[i] = 0.0;
  //   }
  // }
  destroyDouble1DArray(xr);
  destroyDouble1DArray(yr);
  destroyDouble1DArray(zr);
}



void getPolarAngleTrigOverHemisphere(
        int nt,
        double* cosAlpha,
        double* sinAlpha,
        double* cosBeta,
        double* sinBeta,
        double** amp)
{

  int i, j;

  double points = (nt+1) * (2*nt+1);
  double** xr = createDouble2DMatrix( (nt+1), (2*nt+1) );
  double** yr = createDouble2DMatrix( (nt+1), (2*nt+1) );
  double** zr = createDouble2DMatrix( (nt+1), (2*nt+1) );
  
  getDirectionCosineSquareOverHemishpereAndWeights(nt, xr, yr, zr, amp);

  // // OCPolarAngleTrig **trigs = create2DOCPolarAngleTrigArray(nt+1, 2*nt+1);

  // Evaluate sqrt of zr to get cos(beta)
  vdSqrt(points, &zr[0][0], &cosBeta[0]);

   // exaluate A = x + y
  vdAdd(points, &xr[0][0], &yr[0][0], &sinBeta[0]);
  // Take sqrt of A to get sin(beta)
  vdSqrt(points, &sinBeta[0], &sinBeta[0]);

  // Evaluate sqrt of xr
  vdSqrt(points, &xr[0][0], &xr[0][0]);

  // Evaluate sqrt of xr
  vdSqrt(points, &yr[0][0], &yr[0][0]);

  int ii=0;
  for( i = 0; i <= nt; i++) {
    for( j = 0; j <= 2*nt; j++) {
      // cosBeta[i][j] = sqrt(zr[i][j]);
      // sinBeta[i][j] = sqrt(xr[i][j] + yr[i][j]);
      if (sinBeta[ii] != 0.){
        cosAlpha[ii] = xr[i][j] / sinBeta[ii];
        sinAlpha[ii] = yr[i][j] / sinBeta[ii];
        ii++;
      }
      else{
        cosAlpha[ii] = 1.0;
        sinAlpha[ii] = 0.0;
        ii++;
      }
    }
  }
  destroyDouble2DMatrix(xr);
  destroyDouble2DMatrix(yr);
  destroyDouble2DMatrix(zr);
}

// static inline void tent(double freq1, \
//           double freq2, \
//           double freq3, \
//           double amp, \
//           double *spec, \
//           int points, \
//           double fstart, \
//           double finc) {

// // double precision, dimension(0:points-1), intent(inout) :: spec

// double df1, df2, f1, f2, top, t;
// int p, pmid, pmax, i, j;

// double f[3] = {0.0, 0.0, 0.0};
// // double precision, dimension(0:2) :: f


// f[0] = freq1;
// f[1] = freq2;
// f[2] = freq3;

// for( j = 1; j <= 2; j++) {
//     t = f[j];
//     i=j-1;
//     while(i >= 0 && f[i] > t){
//         f[i+1] = f[i];
//         i--;
//     }
//     f[i+1]=t;
// }

// top = amp*2.0 / (f[2]-f[0]);
// p = floor( (f[0]-fstart) / finc );
// pmid = floor( (f[1]-fstart) /finc);
// pmax = floor( (f[2]-fstart) /finc);
// df1 = 0.;
// df2 = 0.;
// if( (f[1]-f[0]) != 0.) df1 = top / (2.0 * (f[1]-f[0]) );
// if( (f[2]-f[1]) != 0.) df2 = top / (2.0 * (f[2]-f[1]) );

// if((pmax < points) && (p >= 0)) {
//     if(p != pmid) {
//         f2 = finc * ( (double)p +1.) + fstart;
//         spec[p] += (f2-f[0]) * (f2-f[0]) * df1;
//         p++;
//         f1 = f2;
//         while(p != pmid) {
//             f2 = finc * ( (double)p + 1.) + fstart;
//             spec[p] += finc * ( (f2-f[0]) + (f1-f[0]) ) * df1;
//             p++;
//             f1 = f2;
//         }
//         spec[p] += (f[1]-f1) * ( (f[1]-f[0]) + (f1-f[0]) ) * df1;
//     } else {
//         spec[p] += (f[1]-f[0]) * top/2.0;
//     }

//     if(p != pmax) {
//         f2 = finc * ( (double)pmid + 1.) + fstart;
//         spec[p] += (f2-f[1]) * ( (f[2]-f2) +(f[2]-f[1]) )* df2;
//         p++;
//         f1 = f2;
//         while(p != pmax) {
//             f2 = finc * ( (double)p + 1.) + fstart;
//             spec[p] += finc * ( (f[2]-f1) + (f[2]-f2) ) * df2;
//             p++;
//             f1 = f2;
//         }
//         spec[p] += ( f[2]-f1 ) * (f[2]-f1) * df2;
//     } else {
//         spec[p] += ( f[2]-f[1] ) * top/2.0;
//     }
//   }
// }



// tent_amp is an optimized version of tent. Using tent2 can provide a
// factor of two boost in the compulation time of spectrum.
// static inline int tent_amp(double *freq1,
//           double *freq2,
//           double *freq3,
//           double *offset,
//           double *amp1,
//           double *amp2,
//           double *amp3,
//           double *spec,
//           int points) {

// double df1, df2, top, t, ampt, diff, Vxy_Vxz, Vyz_Vxz;
// double f10, f21, f20, Vxz, Vxy, Vyz, v1, v2;
// double amp_base;
// int p, pmid, pmax, i, j;

// // off = (int) offset[0];
// p = (int) (freq1[0] + offset[0]);
// if (p == (int)freq2[0] && p == (int)freq3[0]){
//   if(p >= points || p < 0) return 0;
//   spec[p] += (amp1[0]+amp2[0]+amp3[0])/3.0;
// 	return 0;
// }

// double f[3]; //= {0.0, 0.0, 0.0};
// double amp[3];

// f[0] = freq1[0] + offset[0];
// f[1] = freq2[0] + offset[0];
// f[2] = freq3[0] + offset[0];

// amp[0] = amp1[0];
// amp[1] = amp2[0];
// amp[2] = amp3[0];

// for( j = 1; j <= 2; j++) {
//     t = f[j];
//     ampt = amp[j];
//     i=j-1;
//     while(i >= 0 && f[i] > t){
//         f[i+1] = f[i];
//         amp[i+1] = amp[i];
//         i--;
//     }
//     f[i+1]=t;
//     amp[i+1]=ampt;
// }


// p = (int) f[0];
// pmid = (int) f[1];
// pmax = (int) f[2];
// // df1 = 0.;
// // df2 = 0.;
// f10 = f[1]-f[0];
// f21 = f[2]-f[1];
// f20 = f[2]-f[0];

// top = 2.0 / (f[2]-f[0]);

// // if( f10 != 0.) df1 = top / (2.0 * f10 );
// // if( f21 != 0.) df2 = top / (2.0 * f21 );

// // Volume elements pre-factors
// Vxz = (amp[2]-amp[0])/f20;
// Vxy = (amp[1]-amp[0])/f10;
// Vyz = (amp[2]-amp[1])/f21;
// amp_base = amp[0] + Vxz*f10;
// Vxy_Vxz = Vxy+Vxz;
// Vyz_Vxz = Vyz+Vxz;

// if((pmax >= points) || (p < 0)) return 0;

// if(p != pmid) {
//     df1 = top / (2.0 * f10);
//     diff = (double)p +1. - f[0];
//     // cnst = 2.0 * df1 * (amp[0] + diff * (Vxy+Vxz));
//     v1 = diff * diff * df1 * (amp[0] + 0.33333333 * diff * Vxy_Vxz);
//     spec[p++] += v1;
//     diff+=1.0;
//     v2 = diff * diff * df1 * (amp[0] + 0.33333333 * diff * Vxy_Vxz);
//     // vol = v2-v1-cnst;
//     // incr = 2.0 * df1 * (Vxy+Vxz);
//     // i=0;
//     while(p != pmid) {
//       // vol += (cnst+ i*incr);
//       spec[p++] += v2-v1; // vol;
//       // i++;
//       v1=v2;
//       diff+=1.0;
//       v2 = diff * diff * df1 * (amp[0] + 0.33333333 * diff * Vxy_Vxz);
//     }
//     // f1 = (double)p;
//     spec[p] += f10 * top * (amp[0]+amp[1]+amp_base)/6.0 - v1; //vol;
//     // spec[p] += (f[1]-f1) * ( f10 + (f1-f[0]) ) * df1;
// } else {
//     spec[p] += f10 * top * (amp[0]+amp[1]+amp_base)/6.0;
// }

// p=pmax;
// if(p != pmid) {
//     df2 = top / (2.0 * f21 );
//     // f2 = (double)p;
//     diff = f[2] - (double)p;
//     // cnst = 2.0 * df2 * (amp[2] - diff * (Vyz+Vxz));
//     v1 = diff * diff * df2 * (amp[2] - 0.33333333 * diff * Vyz_Vxz);
//     spec[p--] += v1;
//     diff+=1.0;
//     v2 = diff * diff * df2 * (amp[2] - 0.33333333 * diff * Vyz_Vxz);
//     // vol = v2-v1-cnst;
//     // incr = -2.0 * df2 * (Vyz+Vxz);
//     // two_df = 2.0 * df2;
//     // i=0;
//     while(p != pmid) {
//       // vol += (cnst + i*incr);
//       spec[p--] += v2-v1; //vol;
//       // i++;
//       v1=v2;
//       diff+=1.0;
//       v2 = diff * diff * df2 * (amp[2] - 0.33333333 * diff * Vyz_Vxz);
//     }
//     // f1 = (double)p;
//     spec[p] += f21 * top * (amp[2]+amp[1]+amp_base)/6.0 - v1 ; //vol;
// } else {
//     spec[p] += f21 * top * (amp[2]+amp[1]+amp_base)/6.0;
// }
// return 0;
// }



// tent2 is an optimized version of tent. Using tent2 can provide a
// factor of two boost in the computation time of spectrum.
static inline int tent2(double *freq1,
          double *freq2,
          double *freq3,
          double *offset,
          double *amp,
          double *spec,
          int *points
          ) {

double df1, df2, top, t, diff, f10, f21;
int p, pmid, pmax, i, j;

// off = (int) offset[0];
p = (int) (freq1[0] + offset[0]);
if (p == (int)freq2 && p == (int)freq3){
  if(p >= points[0] || p < 0) return 0;
  spec[p] += amp[0];
	return 0;
}

double f[3]; //= {0.0, 0.0, 0.0};

f[0] = freq1[0] + offset[0];
f[1] = freq2[0] + offset[0];
f[2] = freq3[0] + offset[0];

for( j = 1; j <= 2; j++) {
    t = f[j];
    i=j-1;
    while(i >= 0 && f[i] > t){
        f[i+1] = f[i];
        i--;
    }
    f[i+1]=t;
}

top = amp[0]*2.0 / (f[2]-f[0]);
p = (int)f[0];          //floor( f[0] );
pmid = (int)f[1];       //floor( f[1] );
pmax = (int)f[2];       //floor( f[2] );
f10 = f[1]-f[0];
f21 = f[2]-f[1];

if((pmax >= points[0]) || (p < 0)) return 0;

if(p != pmid) {
    df1 = top / f10;
    diff = (double)p + 1. - f[0];
    spec[p++] += 0.5 * diff * diff * df1;
    // diff = (diff - 0.5) * df1;
    diff -= 0.5;
    diff *= df1;
    while(p != pmid) {
      diff += df1;
      spec[p++] += diff;
    }
    spec[p] += (f[1]-(double)p) * ( f10 + ((double)p - f[0]) ) * 0.5 * df1;
} else {
    spec[p] += f10 * top * 0.5;
}

if(p != pmax) {
    df2 = top / f21;
    diff = f[2] - (double)p - 1.;
    spec[p++] += (f21 - diff) * ( diff + f21 ) * 0.5 * df2;
    // diff = (diff + 0.5) * df2;
    diff += 0.5;
    diff *= df2;
    while(p != pmax) {
      diff -= df2;
      spec[p++] += diff;
    }
    spec[p] += pow((f[2] - (double)p), 2) * 0.5 * df2;
} else {
    spec[p] += f21 * top * 0.5;
}
return 0;
};


void powderAverageWithTentingSchemeOverOctant2(
        double *spec,
        double *freq,
        int nt,
        double *amp,
        double *offset,
        int m){

  int i=0, j=0, local_index, n_pts = (nt+1)*(nt+2)/2;
  double amp1, amp2, temp;
  double *amp_address, *freq_address;

  /* Interpolate between frequencies by setting up tents */

    local_index = nt-1;
    amp_address = &amp[nt+1];
    freq_address = &freq[nt+1];

    while(i<n_pts-1){
      temp = amp[i+1] + amp_address[j];
      amp1 = amp[i] + temp;
      
      tent2(&freq[i], &freq[i+1], &freq_address[j], \
                  offset, &amp1, spec, &m);

      if (i<local_index){        
        amp2 = temp + amp_address[j+1];
        tent2(&freq[i+1], &freq_address[j], &freq_address[j+1], \
                  offset, &amp2, spec, &m);
      }
      else{
        local_index=j+nt;
        i++;
      }
      i++; j++;
    }
}

// void powderAverageWithTentingSchemeOverOctant(
//         double *spec,
//         double **powfreq,
//         int nt,
//         double **amp,
//         double *offset,
//         int m){

//   int i, j, local_index;
//   double amp1, amp2, temp;

//   /* Interpolate between frequencies by setting up tents */

//   for (j=0; j<=nt-1; j++){
//     local_index = nt-j-1;
//     for (i=0; i<=local_index; i++){
//       temp = amp[i][j+1] + amp[i+1][j];
//       amp1 = amp[i][j] + temp;
//       tent2(&powfreq[i+1][j], &powfreq[i][j+1], &powfreq[i][j], \
//                   offset, &amp1, spec, &m);

//       if (i<local_index){        
//         amp2 = temp + amp[i+1][j+1];
//         tent2(&powfreq[i+1][j], &powfreq[i][j+1], &powfreq[i+1][j+1], \
//                   offset, &amp2, spec, &m);
//       }
//     }
//   }
// }

void powderAverageWithTentingSchemeOverHemisphere(
        double *spec,
        double **powfreq,
        int nt,
        double **amp,
        double *offset,
        int m){

  int i,j;
  double amp1, amp2, temp;

  /* Interpolate between frequencies by setting up tents */

  for (i=0; i<=nt-1; i++){
    for (j=0; j<=nt-1; j++){
      temp = amp[i][j+1] + amp[i+1][j];
      amp1 = amp[i][j] + temp;
      amp2 = temp + amp[i+1][j+1];

      tent2(&powfreq[i+1][j], &powfreq[i][j+1], &powfreq[i][j], \
                  offset, &amp1, spec, &m);
      tent2(&powfreq[i+1][j], &powfreq[i][j+1], &powfreq[i+1][j+1], \
                  offset, &amp2, spec, &m);


      // tent_amp(&powfreq[i+1][j], &powfreq[i][j+1], &powfreq[i][j], offset, \
      //          &amp[i+1][j], &amp[i][j+1], &amp[i][j], spec, m);
      // tent_amp(&powfreq[i+1][j], &powfreq[i][j+1], &powfreq[i+1][j+1], offset, \
      //          &amp[i+1][j], &amp[i][j+1], &amp[i+1][j+1], spec, m);
    }
  }

  // if (octa == 0){
    for (i=0; i<=nt-1; i++){
      for (j=nt; j<=2*nt-1; j++){
        temp = amp[i][j] + amp[i+1][j+1];
        amp1 = temp + amp[i+1][j];
        amp2 = temp + amp[i][j+1];

        tent2(&powfreq[i][j], &powfreq[i+1][j+1], &powfreq[i+1][j], \
                    offset, &amp1, spec, &m);
        tent2(&powfreq[i][j], &powfreq[i+1][j+1], &powfreq[i][j+1], \
                    offset, &amp2, spec, &m);

        // tent_amp(&powfreq[i][j], &powfreq[i+1][j+1], &powfreq[i][j+1], offset, \
        //          &amp[i][j], &amp[i+1][j+1], &amp[i][j+1], spec, m);
        // tent_amp(&powfreq[i][j], &powfreq[i+1][j+1], &powfreq[i][j+1], offset, \
        //          &amp[i][j], &amp[i+1][j+1], &amp[i][j+1], spec, m);
      }
    }
  // }
};
