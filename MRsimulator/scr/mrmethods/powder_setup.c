
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




// triangle_interpolation is an optimized version of tent. 
int triangle_interpolation(
          double *freq1,
          double *freq2,
          double *freq3,
          double *offset,
          double *amp,
          double *spec,
          int *points
          ) {

double df1, df2, top=0.0, t, diff, f10=0.0, f21=0.0;
int p, pmid, pmax, i, j, clip_right = 0, clip_left = 0;

// off = (int) offset[0];
p = (int) (freq1[0] + offset[0]);
if ((int)freq1[0] == (int)freq2[0] && (int)freq1[0] == (int)freq3[0]){
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

top += (amp[0]*2.0 / (f[2]-f[0]));
p = (int)f[0];          //floor( f[0] );
pmid = (int)f[1];       //floor( f[1] );
pmax = (int)f[2];       //floor( f[2] );
f10 += f[1]-f[0];
f21 += f[2]-f[1];

// if((pmax >= points[0]) || (p < 0)) return 0;

if (pmax < 0) return 0;
if (p > points[0]) return 0;

if(pmax >= points[0]){
  pmax = points[0];
  clip_right = 1;
}

if(pmid >= points[0]){
  pmid = points[0];
  clip_right = 1;
}
    
if (p < 0){
  p = 0;
  clip_left = 1;
}
    
if (pmid < 0){
  pmid = 0;
  clip_left = 1;
}


if(p != pmid) {
    df1 = top / f10;
    diff = (double)p + 1. - f[0];
    if (clip_left == 0){
      spec[p++] += 0.5 * diff * diff * df1;
    }
    else{
      spec[p++] += (diff - 0.5) * df1;
    }
    // diff = (diff - 0.5) * df1;
    diff -= 0.5;
    diff *= df1;
    while(p != pmid) {
      diff += df1;
      spec[p++] += diff;
    }
    if(clip_right == 0){
      spec[p] += (f[1]-(double)p) * ( f10 + ((double)p - f[0]) ) * 0.5 * df1;
    }
} else {
      if(clip_right == 0 && clip_left==0){
        spec[p] += f10 * top * 0.5;
      }
}

if(p != pmax) {
    df2 = top / f21;
    diff = f[2] - (double)p - 1.;

    if (clip_left == 0){
      spec[p++] += (f21 - diff) * ( diff + f21 ) * 0.5 * df2;
    }
    else{
      spec[p++] += (diff + 0.5) * df2;
    }
    // diff = (diff + 0.5) * df2;
    diff += 0.5;
    diff *= df2;
    while(p != pmax) {
      diff -= df2;
      spec[p++] += diff;
    }
    if(clip_right == 0){
      spec[p] += pow((f[2] - (double)p), 2) * 0.5 * df2;
    }
} else {
      if(clip_right == 0){
        spec[p] += f21 * top * 0.5;
      }
  }
  return 0;
}


void powderAverageWithTentingSchemeOverOctant2(
        double *spec,
        double *freq,
        int nt,
        double *amp,
        double *offset,
        int m){

  int i=0, j=0, local_index, n_pts = (nt+1)*(nt+2)/2;
  double amp1=0.0, temp;
  double *amp_address, *freq_address;

  /* Interpolate between frequencies by setting up tents */

    local_index = nt-1;
    amp_address = &amp[nt+1];
    freq_address = &freq[nt+1];

    while(i<n_pts-1){
      temp = amp[i+1] + amp_address[j];
      amp1 = temp;
      amp1 += amp[i];
      
      triangle_interpolation(&freq[i], &freq[i+1], &freq_address[j], offset, &amp1, spec, &m);

      if (i<local_index){  
        amp1 = temp;
        amp1 += amp_address[j+1];
        triangle_interpolation(&freq[i+1], &freq_address[j], &freq_address[j+1], offset, &amp1, spec, &m);
      }
      else{
        local_index=j+nt;
        i++;
      }
      i++; j++;
    }
}


void powderAverageWithTentingSchemeOverHemisphere2(
        double *spec,
        double *freq,
        int nt,
        double *amp,
        double *offset,
        int m){

  int i=0, j=0, local_index, n_pts = (nt+1)*(nt+2)/2;
  double amp1=0.0, temp;
  double *amp_address, *freq_address;

  /* Interpolate between frequencies by setting up tents */
    int l=nt;
    local_index = nt-1;
    amp_address = &amp[4*nt];
    freq_address = &freq[4*nt];

    while(i<n_pts-1){
      temp = amp[i+1] + amp_address[j];
      amp1 = temp;
      amp1 += amp[i];
      
      triangle_interpolation(&freq[i], &freq[i+1], &freq_address[j], offset, &amp1, spec, &m);

      if (i<local_index){  
        amp1 = temp;
        if(j==4*l){
          amp1 += amp_address[j-4*l];
          triangle_interpolation(&freq[i+1], &freq_address[j], &freq_address[0], offset, &amp1, spec, &m);
        }
        else amp1 += amp_address[j+1];
        triangle_interpolation(&freq[i+1], &freq_address[j], &freq_address[j+1], offset, &amp1, spec, &m);
        i++; j++;
      }
      else{
        local_index=i+nt;
        i++;
      }
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
//       triangle_interpolation(&powfreq[i+1][j], &powfreq[i][j+1], &powfreq[i][j], \
//                   offset, &amp1, spec, &m);

//       if (i<local_index){        
//         amp2 = temp + amp[i+1][j+1];
//         triangle_interpolation(&powfreq[i+1][j], &powfreq[i][j+1], &powfreq[i+1][j+1], \
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

      triangle_interpolation(&powfreq[i+1][j], &powfreq[i][j+1], &powfreq[i][j], \
                  offset, &amp1, spec, &m);
      triangle_interpolation(&powfreq[i+1][j], &powfreq[i][j+1], &powfreq[i+1][j+1], \
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

        triangle_interpolation(&powfreq[i][j], &powfreq[i+1][j+1], &powfreq[i+1][j], \
                    offset, &amp1, spec, &m);
        triangle_interpolation(&powfreq[i][j], &powfreq[i+1][j+1], &powfreq[i][j+1], \
                    offset, &amp2, spec, &m);

        // tent_amp(&powfreq[i][j], &powfreq[i+1][j+1], &powfreq[i][j+1], offset, \
        //          &amp[i][j], &amp[i+1][j+1], &amp[i][j+1], spec, m);
        // tent_amp(&powfreq[i][j], &powfreq[i+1][j+1], &powfreq[i][j+1], offset, \
        //          &amp[i][j], &amp[i+1][j+1], &amp[i][j+1], spec, m);
      }
    }
  // }
};

// static inline double two_times_triangle_area(
//                         double *a, 
//                         double *b,
//                         double *c)
// {
//     return ((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]));
// }


void rasterization(double * grid,
                   double *v0,
                   double *v1,
                   double *v2,
                   int rows,
                   int columns){


  double A12, B12, C12, A20, B20, C20, A01, B01, C01;
  double minX, minY, maxX, maxY, w0, w1, w2;
  double w0_row, w1_row, w2_row;
  int i, j, i_, minX_, minY_, maxX_, maxY_;

  minX = fmin(fmin(v0[0], v1[0]), v2[0]);
  minY = fmin(fmin(v0[1], v1[1]), v2[1]);
  maxX = fmax(fmax(v0[0], v1[0]), v2[0]);
  maxY = fmax(fmax(v0[1], v1[1]), v2[1]);


  // clip against screen bounds
  minX_ = (int) fmax(minX, 0.);
  minY_ = (int) fmax(minY, 0.);
  maxX_ = (int) fmin(maxX, (double) rows - 1.);
  maxY_ = (int) fmin(maxY, (double) columns - 1.);

  
  A12 = (v2[0]-v1[0]); B12 = (v2[1]-v1[1]); C12 = -A12*v1[1] + B12*v1[0];
  A20 = (v0[0]-v2[0]); B20 = (v0[1]-v2[1]); C20 = -A20*v2[1] + B20*v2[0];
  A01 = (v1[0]-v0[0]); B01 = (v1[1]-v0[1]); C01 = -A01*v0[1] + B01*v0[0];
      

  w0_row = A12 *minY - B12*minX + C12;
  w1_row = A20 *minY - B20*minX + C20;
  w2_row = A01 *minY - B01*minX + C01;
    
  // Rasterize
  for(i=minY_; i<=maxY_; i++){

    // Determine barycentric coordinates
    w0 = w0_row;
    w1 = w1_row;
    w2 = w2_row;

    i_ = rows*i;
    for(j=minX_; j<=maxX_; j++){
      // If p is on or inside all edges, render pixel.
      if((int) w0>=0 && (int) w1>=0 && (int) w2>= 0){
        grid[i_+j] += 1.; //(w0+w1+w2);
      }
      if ((int) w0<=0 && (int) w1<=0 && (int) w2<=0){
        grid[i_+j] += -1.; //(w0+w1+w2);
      }
      // i_++;
              
      w0 -= B12;
      w1 -= B20;
      w2 -= B01;

    }
    w0_row += A12;
    w1_row += A20;
    w2_row += A01;
  }
}

// static inline void get_KL(
//             double *v0,
//             double *v1,
//             double *x,
//             double *y)
// {
//   double K[2] = 
// }
// def get_KL(self, x0, y0, x1, y1, eps = 1e-5):
//     self.x0, self.y0, self.x1, self.y1 = x0, y0, x1, y1
//     Kx, Ky, Lx, Ly = 0, 0, 0, 0
//     for sec in self.clip(0, 1, 1, 0):
//         v1 = Point(sec[0], sec[1])
//         v0 = Point(sec[2], sec[3])
//         if abs(v0.x - 1) < eps and abs(v1.x - 1) < eps \
//         or abs(v0.y - 1) < eps and abs(v1.y - 1) < eps:
//             continue

//         Kx += 1./4 * (v0.y - v1.y)
//         Ky += 1./4 * (v1.x - v0.x)
//         Lx += 1./8 * (v0.y-v1.y) * (v0.x+v1.x)
//         Ly += 1./8 * (v1.x-v0.x) * (v0.y+v1.y)
//     return Point(Kx, Ky), Point(Lx, Ly)


// static inline void clip(
//             double left,
//             double right,
//             double bottom,
//             double top,
//             double *p0,
//             double *p1,
//             double *clipped_lines
//             )
// {
//   int edge, i;
//   double t0=0, t1=1, r;
//   double delta_x = p1[0] - p0[0], delta_y=p1[1]-p0[1], p, q;

//   for(edge=0; edge<4; edge++){
//     if(edge ==0){
//       p = -delta_x;
//       q = -(left - p0[0]);
//     }
//     else if(edge ==1){
//       p = delta_x;
//       q = (right - p0[0]);
//     }
//     else if(edge ==2){
//       p = delta_y;
//       q = (bottom - p0[1]);
//     }
//     else if(edge ==3){
//       p = -delta_y;
//       q = -(top - p0[1]);
//     }
//     if(p == 0 && q < 0) return 0;
//     if (p < 0){
//       r = q / (double)p;
//       if(r > t1)  return 0;
//       if(r > t0) t0 = r;  // line is clipped!
//     }
//     else if(p > 0){
//       r = q / (double)p;
//       if(r < t0) return 0;
//       if(r < t1)  t1 = r;   // line is clipped!
//     }
//   }
// }

//     def clip(self, left, right, bottom, top):
//         t0, t1 = 0, 1
//         xdelta = self.x1 - self.x0
//         ydelta = self.y1 - self.y0
//         for edge in range(4): #traverse through left, right, bottom, top edges.
//             if   edge == 0:   p, q = -xdelta, -(left-self.x0) 
//             elif edge == 1:   p, q =  xdelta,  (right-self.x0)
//             elif edge == 2:   p, q =  ydelta,  (bottom-self.y0)
//             elif edge == 3:   p, q = -ydelta, -(top-self.y0)
//             if p == 0 and q < 0:    return []
//             if p < 0:
//                 r = q / float(p)
//                 if r > t1:          return []
//                 elif r > t0:        t0 = r   # line is clipped!
//             elif p > 0:
//                 r = q / float(p)
//                 if r < t0:          return []
//                 elif r < t1:        t1 = r   # line is clipped!
//         clipped_line = (self.x0 + t0*xdelta, self.y0 + t0*ydelta, 
//                         self.x0 + t1*xdelta, self.y0 + t1*ydelta)
//         return [clipped_line]