/* basic Monte Carlo simulation */
#include "lj_npr.h"
#include "av.h"
#include "bpav.h"



int n = 108;
int nequil = 100000;
int nsteps = 1000000;
double rho = 0.62;
double tp = 1.15;
double rcdef = 1e9; /* half-box cutoff */
const char *fnpos = "lj.pos";



/* energy increase of a new position */
double lj_widom(lj_t *lj, double beta)
{
  double xi[D], dx[D], l = lj->l, invl = 1/l;
  double dr2, ir2, ir6, de, ep, ep6 = 0, ep12 = 0, eptail = 0;
  int d, j, n = lj->n;

  /* randomly choose a position */
  for ( d = 0; d < D; d++ ) xi[d] = rand01() * lj->l;

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
  de = ep + eptail - lj->epot_tail;

  return exp(-beta * de);
}



int main(void)
{
  int t, acc;
  lj_t *lj;
  double epsm = 0, psm = 0, accsm = 0, widsm = 0;
  double mu, beta = 1/tp, bp, dbp, amp = 0.1/rho;
  bpav_t bpav[1];

  mtscramble(time(NULL));
  bpav_clear(bpav);
  lj = lj_open(n, rho, rcdef);
  lj->dof = n * D;
  //lj->vol -= 0.5; lj_setrho(lj, n / lj->vol); 
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
    widsm += lj_widom(lj, beta);
    bpav_add(bpav, lj, beta);
  }
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
  mu = tp * log((n+1)/lj->vol) - tp * log(widsm/nsteps);
  bp = bpav_get(bpav, &dbp);
  printf("rho %g, tp %g, ep %g, p %g, mu %g, acc %g%%\n",
      rho, tp, epsm/nsteps/n, psm/nsteps,
      mu, 100.*accsm/nsteps);
  printf("vol %g, bp %g, dbp %g\n", lj->vol, bp, dbp);
  return 0;
}

