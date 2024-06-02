#ifndef __tables__
#define __tables__

#define lerp_plus(x, i) lerp((w), gauss_table[(i)], gauss_table[(i) + 1])
#define lerp_minus(x, i) lerp((w), gauss_table[(i)], gauss_table[(i) - 1])

/** Cosine and Sine tables **/
double cos_table[62833];
double sin_table[62833];
double trig_table_precision_inverse;

static inline void generate_trig_table(void) {
  extern double cos_table[62833];
  extern double sin_table[62833];
  extern double trig_table_precision_inverse;

  int i, n = 62833;
  double precision = 0.0001;
  trig_table_precision_inverse = 1.0 / precision;

  for (i = 0; i < n - 1; i++) {
    cos_table[i] = cos(i * precision);
    sin_table[i] = sin(i * precision);
  }

  cos_table[n - 1] = 1.0;
  sin_table[n - 1] = 0.0;
}

/** Gaussian tables **/
double gauss_table[1500];
double gauss_table_precision_inverse;

static inline void generate_gauss_table(void) {
  extern double gauss_table[1500];
  extern double gauss_table_precision_inverse;

  int i, n = 1500;
  double precision = 0.002, sigma = 1.0 / 4.0, factor, temp;
  gauss_table_precision_inverse = 1.0 / precision;

  factor = 1.0 / (2.0 * sigma * sigma);
  for (i = 0; i <= n - 1; i++) {
    temp = (i * precision);
    gauss_table[i] = exp(-temp * temp * factor);
  }
}

static inline void generate_tables(void) {
  generate_trig_table();
  generate_gauss_table();
}

#endif /* __tables__ */
