/* Monte Carlo simulation with harmonic restraints */
#include <time.h>
#include "av.h"
#include "lj.h"



int nmax = 140; // maximal number of particles
double nr = 107.5; // average number
double knr = 10.0; // spring constant of the number
double vol = 150; // average volume
double kvol = 10.0;
int nequil = 100000;
int nsteps = 1000000;
double rho = 0.9;
double tp = 1.15;
double rcdef = 1e9; /* half-box cutoff */
double amp = 0.2; /* Monte Carlo move size */
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



/* volume move */
int lj_vmove(lj_t *lj, double beta, double vol, double kvol)
{
  const double lnvamp = 0.01;
  double lnlo, lnln, lo, ln, vo, vn, s, ss, is6;
  double epo, epn, ep6n, ep12n, eptailn, dw;
  int i, n = lj->n, acc;

  vo = lj->vol;
  lo = lj->l;
  lnlo = log(lo);
  lnln = lnlo + lnvamp/D * (rand01() * 2 - 1);
  ln = exp( lnln );
  vn = exp( D * lnln );

  epo = lj->epot;
  s = ln/lo;
  ss = s * s;
  is6 = 1 / (ss * ss * ss);
  ep6n = lj->ep6 * is6;
  ep12n = lj->ep12 * is6 * is6;
  /* assume half-box cutoff here */
  eptailn = lj_gettail(ln/2, lj->n / vn, lj->n, NULL);
  epn = ep12n - ep6n + eptailn;

  dw = beta * (epn - epo) + (lnln - lnlo) * D * n +
     - kvol * (vn - vo) * (vn + vo - 2 * vol);
  acc = 1;
  if ( dw < 0 ) {
    double r = rand01();
    acc = ( r < exp(dw) );
  }
  if ( acc ) {
    lj_setrho(lj, lj->n / vn);
    lj->ep12 = ep12n;
    lj->ep6 = ep6n;
    lj->ep0 = ep12n - ep6n;
    lj->vir = ep12n * 12 - ep6n * 6;
    lj->epot = lj->ep0 + lj->epot_tail;
    /* scale the coordinates */
    for ( i = 0; i < lj->n; i++ ) {
      vsmul(lj->x[i], s);
    }
    for ( i = 0; i < lj->n * lj->n; i++ ) {
      lj->r2ij[i] *= ss;
    }
  }
  return acc;
}


/* particle move */
int lj_nmove(lj_t *lj, double beta, double nr, double knr)
{
  double xi[D], dx[D], l = lj->l, invl = 1/l, dw;
  double dr2, ir2, ir6, de, ep, ep6 = 0, ep12 = 0, eptail = 0;
  int d, i, j, n = lj->n, acc;

  if ( rand01() < 0.5 ) { /* try to add a particle */
    if ( n + 1 > nmax ) return 0;

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

    dw = log( lj->vol / (n + 1) ) - beta * de
       + knr * (2 * nr - 2 * n - 1);
    acc = ( dw > 0 || rand01() < exp(dw) );
    if ( acc ) {
      vcopy(lj->x[n], xi);
      // ignore v and f
      lj->n = n + 1;
      lj_setrho(lj, lj->n/lj->vol);
      lj->dof = lj->n * D;
      //printf("de %g\n", de);
    }

  } else { /* try to remove a particle */
    if ( n == 0 ) return 0;

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

    dw = log( n / lj->vol ) + beta * de
       - knr * (2 * nr - 2 * n + 1);
    acc = ( dw > 0 || rand01() < exp(dw) ) * (-1);
    if ( acc ) {
      /* removing particle i */
      if ( i < n - 1 ) {
        vcopy(lj->x[i], lj->x[n - 1]);
        // ignore v and f
      }
      lj->n = n - 1;
      lj_setrho(lj, lj->n/lj->vol);
      lj->dof = lj->n * D;
    }
  }

  return acc;
}



int main(void)
{
  int t, acc, nacc = 0, vacc = 0;
  lj_t *lj;
  av_t avep[1], avp[1], avacc[1], avnr[1], avvol[1], avwid[1];
  double beta = 1/tp, bmu = 0;

  mtscramble(time(NULL));
  av_clear(avep);
  av_clear(avp);
  av_clear(avacc);
  av_clear(avnr);
  av_clear(avvol);
  av_clear(avwid);
  lj = lj_open(nmax, rho, rcdef);
  lj_energy(lj);
  /* equilibration */
  for ( t = 1; t <= nequil; t++ )
    lj_metro(lj, amp, beta);
  /* production */
  for ( t = 1; t <= nsteps; t++ ) {
    acc = lj_metro(lj, amp, beta);
    if ( t % 10 == 0 ) {
      nacc = lj_nmove(lj, beta, nr, knr);
      av_add(avnr, lj->n);
      vacc = lj_vmove(lj, beta, vol, kvol);
      av_add(avvol, lj->vol);
      if ( t % 100000 == 0 ) {
        printf("t %d, n %d, nr %g, %g, nacc %d, "
               "vol %g, av. vol %g, %g, vacc %d\n",
           t, lj->n, av_getave(avnr), nr, nacc,
           lj->vol, av_getave(avvol), vol, vacc);
      }
    }
    av_add(avep, lj->epot);
    av_add(avp, lj_calcp(lj, tp));
    av_add(avacc, acc);
    av_add(avwid, lj_widom(lj, beta));
  }
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
  printf("tp %g, ep %g, nr %g, vol %g, acc %g%%, p %g, bmu %g\n",
      tp, av_getave(avep)/nr, av_getave(avp),
      av_getave(avnr), av_getave(avvol),
      100.*av_getave(avacc), av_getave(avp),
      -log(av_getave(avwid)*lj->vol/nr) );
  return 0;
}

