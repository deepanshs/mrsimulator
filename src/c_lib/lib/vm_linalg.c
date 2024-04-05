
#include "vm_linalg.h"

// Find the roots of cubic equation
void vm_haeberlen_components(int n, double *expr_base_p, double *expr_base_q,
                             double zeta, double eta, double rho, double *param) {
  int counter = n;
  bool status;
  double z2, ze, z2e2, z3e2, z2e, r2, p_prime, q_prime;
  double root_0, root_2, root_1, temp, arg, a_cos, angle;

  angle = 2.09439510239;

  z2 = zeta * zeta;
  ze = zeta * eta;
  z2e2 = ze * ze;
  z2e = z2 * eta;
  z3e2 = zeta * z2e2;
  r2 = rho * rho;

  double basis_rho_q[8] = {rho * r2, zeta * r2,  z2 * rho, z2 * zeta,
                           ze * r2,  z2e2 * rho, z3e2,     z2e * rho};
  double basis_rho_p[5] = {r2, zeta * rho, z2, ze * rho, z2e2};

  while (counter-- > 0) {
    p_prime = cblas_ddot(5, basis_rho_p, 1, expr_base_p++, n);
    q_prime = cblas_ddot(8, basis_rho_q, 1, expr_base_q++, n);

    temp = sqrt(p_prime);
    arg = q_prime / (p_prime * temp);
    a_cos = acos(arg) / 3.0;
    status = a_cos < 0.5235987756;

    root_0 = get_cos_from_table(a_cos);
    a_cos -= angle;
    root_1 = get_cos_from_table(a_cos);
    a_cos -= angle;
    root_2 = get_cos_from_table(a_cos);

    temp *= 2.0;
    if (status) {
      *param++ = temp * root_0;
      *param++ = (root_1 - root_2) / root_0;
    } else {
      *param++ = temp * root_2;
      *param++ = (root_1 - root_0) / root_2;
    }
  }
}
