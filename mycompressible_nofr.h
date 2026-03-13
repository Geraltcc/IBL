/**
# Compressible gas dynamics (no flux redistribution version)
*/

#include "myconservation_nofr.h"

scalar rho[], E[];
vector w[];
scalar * scalars = {rho, E};
vector * myvectors = {w};
double gammao = 1.4;

void flux (const double * s, double * f, double * e)
{
  double rho = s[0], E = s[1], wn = s[2], w2 = 0.;
  for (int i = 2; i < 2 + dimension; i++)
    w2 += sq(s[i]);

  if (rho <= 0. || rho != rho) {
    for (int i = 0; i < 2 + dimension; i++)
      f[i] = 0.;
    e[0] = e[1] = 0.;
    return;
  }

  double un = wn/rho, p = (gammao - 1.)*(E - 0.5*w2/rho);

  if (p <= 0. || p != p) {
    f[0] = wn;
    f[1] = un*E;
    f[2] = un*wn;
    for (int i = 3; i < 2 + dimension; i++)
      f[i] = un*s[i];
    e[0] = un;
    e[1] = un;
    return;
  }

  f[0] = wn;
  f[1] = un*(E + p);
  f[2] = un*wn + p;
  for (int i = 3; i < 2 + dimension; i++)
    f[i] = un*s[i];

  double c = sqrt(gammao*p/rho);
  e[0] = un - c;
  e[1] = un + c;
}

event cell_linking(i++) {
  for (scalar s in evolving) {
    foreach() {
      if (cs[] > 0. && cs[] < CS_FLOOR) {
        double sum_val = 0., sum_w = 0.;
        foreach_neighbor(1) {
          if (cs[] >= CS_FLOOR) {
            sum_val += cs[] * s[];
            sum_w += cs[];
          }
        }
        if (sum_w > 0.)
          s[] = sum_val / sum_w;
      }
    }
  }
}