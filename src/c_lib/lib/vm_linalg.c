
#include "vm_linalg.h"

static void inline pas_to_haeberlen(double r1, double r2, double r3, double *param) {
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
void vm_haeberlen_components(int n, double *expr_base_p, double *expr_base_q,
                             double zeta, double eta, double rho, double *param) {
  int counter = n;
  double z2, ze, z2e2, z3, z3e2, z2e, r2, p, q;
  double root_1, root_2, root_3, temp, arg, a_cos, angle;

  angle = CONST_PI * 0.6666666666666666;

  z2 = zeta * zeta;
  ze = zeta * eta;
  z3 = z2 * zeta;
  z2e2 = ze * ze;
  z2e = z2 * eta;
  z3e2 = zeta * z2e2;
  r2 = rho * rho;

  double basis_rho_q[8] = {rho * r2, zeta * r2,  z2 * rho, z3,
                           ze * r2,  z2e2 * rho, z3e2,     z2e * rho};
  double basis_rho_p[8] = {r2, zeta * rho, z2, z3, ze * rho, z2e2, z3e2, z2e};

  while (counter-- > 0) {
    p = cblas_ddot(8, basis_rho_p, 1, expr_base_p++, n);
    q = cblas_ddot(8, basis_rho_q, 1, expr_base_q++, n);

    temp = sqrt(-1.0 / p) * 0.8660254037844386;  // np.sqrt(3) / 2
    arg = (temp * 3.0 * q) / p;
    a_cos = acos(arg) / 3.0;

    root_1 = get_cos_from_table(a_cos) / temp;
    root_2 = get_cos_from_table(a_cos - angle) / temp;
    root_3 = get_cos_from_table(a_cos + angle) / temp;

    pas_to_haeberlen(root_1, root_2, root_3, param);
    param += 2;
  }
}
