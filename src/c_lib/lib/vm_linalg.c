
#include "vm_linalg.h"

// Find the roots of cubic equation
void vm_haeberlen_components(int n, double *expr_base_p, double *expr_base_q,
                             double zeta, double eta, double rho, double *param) {
  int counter = n;
  bool status;
  double z2, ze, e2, z2e, r2, p_prime, q_prime;
  double root_0, root_2, root_1, temp, arg, a_cos, angle;

  angle = 2.09439510239;

  z2 = zeta * zeta;
  e2 = eta * eta;
  ze = zeta * eta;
  z2e = z2 * eta;
  r2 = rho * rho;

  double basis_rho_q[6] = {
      rho * r2, zeta * r2, ze * r2, z2 * rho * (-3.0 + e2), z2 * zeta * (1.0 - e2),
      z2e * rho};
  double basis_rho_p[4] = {r2, zeta * rho, ze * rho, z2 * (3.0 + e2)};

  while (counter-- > 0) {
    p_prime = cblas_ddot(4, basis_rho_p, 1, expr_base_p, n);
    q_prime = cblas_ddot(6, basis_rho_q, 1, expr_base_q, n);

    expr_base_p++;
    expr_base_q++;

    temp = sqrt(p_prime);
    arg = q_prime / (p_prime * temp);
    a_cos = acos(arg) / 3.0;
    status = a_cos < 0.5235987756;

    root_0 = cos(a_cos);
    a_cos -= angle;
    root_1 = cos(a_cos);
    a_cos -= angle;
    root_2 = cos(a_cos);

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