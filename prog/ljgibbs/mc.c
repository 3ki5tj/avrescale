/* basic Monte Carlo simulation */
#include "lj.h"



int n = 108;
int nequil = 100000;
int nsteps = 1000000;
double rho = 0.08;
double tp = 1.15;
double rcdef = 1e9; /* half-box cutoff */
double amp = 0.2; /* Monte Carlo move size */
const char *fnpos = "lj.pos";



/* energy increase of a new position */
double lj_denp(lj_t *lj, double *xi)
{
  double dx[D], l = lj->l, invl = 1/l;
  double dr2, ir2, ir6, ep, ep6 = 0, ep12 = 0, eptail = 0;
  int j, n = lj->n;

  /* compute the energy change */
  for ( j = 0; j < n; j++ ) {
    dr2 = lj_pbcdist2(dx, xi, lj->x[j], l, invl);
    if ( dr2 >= lj->rc2 ) continue;
    ir2 = 1 / dr2;
    ir6 = ir2 * ir2 * ir2;
    ep12 += ir6 * ir6;
    ep6 += ir6;
  }
  ep6 *= 4;
  ep12 *= 4;
  ep = ep12 - ep6;

  /* compute the change of tail correction */
  eptail = lj_gettail(lj->rc, (n + 1)/lj->vol, n, NULL);
  return ep + eptail - lj->epot_tail;
}

double lj_widom(lj_t *lj, double beta)
{
  double xi[D], de;
  int d;

  /* randomly choose a position */
  for ( d = 0; d < D; d++ ) xi[d] = rand01() * lj->l;

  de = lj_denp(lj, xi);
  return exp(-beta*de);
}

int main(void)
{
  int t, acc;
  lj_t *lj;
  double epsm = 0, psm = 0, accsm = 0, widsm = 0;
  double mu;

  mtscramble(time(NULL));
  lj = lj_open(n, rho, rcdef);
  lj_energy(lj);
  /* equilibration */
  for ( t = 1; t <= nequil; t++ )
    lj_metro(lj, amp, 1/tp);
  /* production */
  for ( t = 1; t <= nsteps; t++ ) {
    acc = lj_metro(lj, amp, 1/tp);
    epsm += lj->epot;
    psm += lj_calcp(lj, tp);
    accsm += acc;
    widsm += lj_widom(lj, 1/tp);
  }
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
  mu = tp * log((n+1)/lj->vol) - tp * log(widsm/nsteps);
  printf("rho %g, tp %g, ep %g, p %g, mu %g, acc %g%%\n",
      rho, tp, epsm/nsteps/n, psm/nsteps,
      mu, 100.*accsm/nsteps);
  return 0;
}

