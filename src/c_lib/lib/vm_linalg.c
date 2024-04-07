
#include "vm_linalg.h"

// Find the roots of cubic equation and then zeta, eta based on Haeberlen convention
void vm_haeberlen_components(int n, double *expr_base_p, double *expr_base_q,
                             double zeta, double eta, double rho, double *param) {
  int counter = n, index;
  bool status;
  double *param_start = param;
  double z2, zer, e2, r2, p_prime, q_prime;
  double root_0, root_2, root_1, temp, arg, a_cos, angle;

  angle = 2.09439510239 * trig_table_precision_inverse;  // 2 pi / 3

  z2 = zeta * zeta;
  e2 = eta * eta;
  zer = zeta * eta * rho;
  r2 = rho * rho;

  double basis_rho_q[6] = {
      rho * r2,  zeta * r2, zer * rho, z2 * rho * (-3.0 + e2), z2 * zeta * (1.0 - e2),
      zeta * zer};
  double basis_rho_p[4] = {r2, zeta * rho, zer, z2 * (3.0 + e2)};

  while (counter-- > 0) {
    p_prime = cblas_ddot(4, basis_rho_p, 1, expr_base_p++, n);
    q_prime = cblas_ddot(6, basis_rho_q, 1, expr_base_q++, n);

    temp = sqrt(p_prime);
    arg = q_prime / (p_prime * temp);
    a_cos = acos(arg) / 3.0;
    status = a_cos < 0.5235987756;  // pi / 6

    a_cos *= trig_table_precision_inverse;
    index = (int)a_cos;
    root_0 = lerp(a_cos - index, cos_table[index], cos_table[index + 1]);

    a_cos += angle;
    index = (int)a_cos;
    root_2 = lerp(a_cos - index, cos_table[index], cos_table[index + 1]);

    a_cos += angle;
    index = (int)a_cos;
    root_1 = lerp(a_cos - index, cos_table[index], cos_table[index + 1]);

    if (status) {
      *param++ = temp * root_0;
      *param++ = (root_1 - root_2) / root_0;
    } else {
      *param++ = temp * root_2;
      *param++ = (root_1 - root_0) / root_2;
    }
  }
  cblas_dscal(n, 2.0, param_start, 2);
}
