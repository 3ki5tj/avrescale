/* Monte Carlo simulation with fractional number of particles */
#include <time.h>
#include "av.h"
#include "lj.h"



double nr = 107.5;
int nequil = 100000;
int nsteps = 1000000;
double rho = 0.08;
double tp = 1.15;
double rcdef = 1e9; /* half-box cutoff */
double amp = 0.2; /* Monte Carlo move size */
const char *fnpos = "lj.pos";
int adaptive = 1;



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



/* particle move */
int lj_nmove(lj_t *lj, double beta, double bmu)
{
  double xi[D], dx[D], l = lj->l, invl = 1/l, dw;
  double dr2, ir2, ir6, de, ep, ep6 = 0, ep12 = 0, eptail = 0;
  int d, i, j, n = lj->n, acc;

  if ( nr > n ) { /* try to add a particle */
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

    dw = log( lj->vol / (n + 1) ) + bmu - beta * de;
    acc = ( dw > 0 || rand01() < exp(dw) );
    if ( acc ) {
      vcopy(lj->x[n], xi);
      // ignore v and f
      lj->n = n + 1;
      lj_setrho(lj, lj->n/lj->vol);
      //printf("de %g\n", de);
    }

  } else { /* try to remove a particle */
    /* randomly choose a particle to remove */
    i = (int) (rand01() * n);

    /* compute the energy change */
    for ( j = 0; j < n; j++ ) {
      if ( j == i ) continue;
      dr2 = lj_pbcdist2(dx, lj->x[i], lj->x[j], l, invl);
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
    eptail = lj_gettail(lj->rc, (n - 1)/lj->vol, n, NULL);
    de = ep + lj->epot_tail - eptail;

    dw = log( n / lj->vol ) - bmu + beta * de;
    acc = ( dw > 0 || rand01() < exp(dw) );
    if ( acc ) {
      /* removing particle i */
      if ( i < n - 1 ) {
        vcopy(lj->x[i], lj->x[n - 1]);
        // ignore v and f
      }
      lj->n = n - 1;
      lj_setrho(lj, lj->n/lj->vol);
    }
  }

  return acc;
}



int main(void)
{
  int t, acc;
  lj_t *lj;
  av_t avep[1], avp[1], avacc[1], avnr[1], avwid[1];
  double beta = 1/tp, bmu = 0;

  mtscramble(time(NULL));
  av_clear(avep);
  av_clear(avp);
  av_clear(avacc);
  av_clear(avnr);
  av_clear(avwid);
  lj = lj_open((int) nr + 1, rho, rcdef);
  lj_energy(lj);
  /* equilibration */
  for ( t = 1; t <= nequil; t++ )
    lj_metro(lj, amp, beta);
  /* production */
  for ( t = 1; t <= nsteps; t++ ) {
    acc = lj_metro(lj, amp, beta);
    if ( t % 10 == 0 ) {
      lj_nmove(lj, beta, bmu);
      av_add(avnr, lj->n);
      if ( adaptive ) {
        bmu += -0.001*(lj->n - nr);
        if ( t % 100 == 0 ) {
          printf("t %d, bmu %g, n %d, nr %g, %g\n",
             t, bmu, lj->n, av_getave(avnr), nr);
        }
      }
    }
    av_add(avep, lj->epot);
    av_add(avp, lj_calcp(lj, tp));
    av_add(avacc, acc);
    av_add(avwid, lj_widom(lj, beta));
  }
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
  printf("rho %g, tp %g, ep %g, p %g, nr %g, acc %g%%, bmu %g\n",
      rho, tp, av_getave(avep)/nr, av_getave(avp),
      av_getave(avnr), 100.*av_getave(avacc),
      -log(av_getave(avwid)*lj->vol/nr) );
  return 0;
}

