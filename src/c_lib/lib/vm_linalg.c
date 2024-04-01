
#include "vm_linalg.h"

static void inline zeta_eta(double r1, double r2, double r3, double *param) {
  double ar1, ar2, ar3;

  ar1 = absd(r1);
  ar2 = absd(r2);
  ar3 = absd(r3);

  if (ar1 > ar2 && ar1 > ar3) {
    *param++ = r1;                  //  zeta
    *param = absd((r2 - r3) / r1);  //  eta
  } else if (ar2 > ar1 && ar2 > ar3) {
    *param++ = r2;                  //  zeta
    *param = absd((r1 - r3) / r2);  //  eta
  } else {
    *param++ = r3;                  //  zeta
    *param = absd((r1 - r2) / r3);  //  eta
  }
}

// Find the roots of cubic equation
void cubic_roots(int n, double *expr_base_p, double *expr_base_q, double zeta,
                 double eta, double rho, double *param) {
  int counter = n;
  double z2, ze, z2e2, z3, z3e2, z2e, r2, r3, p, q;
  double root_1, root_2, root_3, temp, arg, a_cos, angle;
  double *basis_rho_q, *basis_rho_p;
  basis_rho_q = malloc_double(8);
  basis_rho_p = malloc_double(8);

  angle = CONST_PI * 0.6666666666666666;

  while (counter-- > 0) {
    z2 = zeta * zeta;
    ze = zeta * eta;
    z3 = z2 * zeta;
    z2e2 = ze * ze;
    z2e = z2 * eta;
    z3e2 = zeta * z2e2;
    r2 = rho * rho;
    r3 = rho * r2;

    // [3, 2, 1, 0, 2, 1, 0, 1]
    basis_rho_q[0] = r3;
    basis_rho_q[1] = zeta * r2;
    basis_rho_q[2] = z2 * rho;
    basis_rho_q[3] = z3;
    basis_rho_q[4] = ze * r2;
    basis_rho_q[5] = z2e2 * rho;
    basis_rho_q[6] = z3e2;
    basis_rho_q[7] = z2e * rho;

    // [2, 1, 0, 0, 1, 0, 0, 0]
    basis_rho_p[0] = r2;
    basis_rho_p[1] = zeta * rho;
    basis_rho_p[2] = z2;
    basis_rho_p[3] = z3;
    basis_rho_p[4] = ze * rho;
    basis_rho_p[5] = z2e2;
    basis_rho_p[6] = z3e2;
    basis_rho_p[7] = z2e;

    p = cblas_ddot(8, basis_rho_p, 1, expr_base_p, n);
    q = cblas_ddot(8, basis_rho_q, 1, expr_base_q, n);

    temp = sqrt(-1.0 / p) * 0.8660254037844386;  // np.sqrt(3) / 2
    arg = (temp * 3.0 * q) / p;
    a_cos = acos(arg) / 3.0;

    root_1 = get_cos_from_table(a_cos) / temp;
    root_2 = get_cos_from_table(a_cos - angle) / temp;
    root_3 = get_cos_from_table(a_cos + angle) / temp;

    zeta_eta(root_1, root_2, root_3, param);
    param += 2;
    expr_base_p++;
    expr_base_q++;
  }
  free(basis_rho_q);
  free(basis_rho_p);
}
