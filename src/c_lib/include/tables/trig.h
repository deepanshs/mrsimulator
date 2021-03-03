#ifndef __tables__
#define __tables__

double cos_table[62833];
double sin_table[62833];
double table_precision_inverse;

static inline void generate_table() {
  extern double cos_table[62833];
  extern double sin_table[62833];
  extern double table_precision_inverse;

  int i, n = 62833;
  double precision = 0.0001;
  table_precision_inverse = 1.0 / precision;

  for (i = 0; i < n - 1; i++) {
    cos_table[i] = cos(i * precision);
    sin_table[i] = sin(i * precision);
  }

  cos_table[n - 1] = 1.0;
  sin_table[n - 1] = 0.0;
}

#endif /* __tables__ */
