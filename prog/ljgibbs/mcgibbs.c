/* Monte Carlo simulation with harmonic restraints */
#include <time.h>
#include "av.h"
#include "lj.h"



int nmax = 108; // maximal number of particles
double cnr[2] = {108, 108}; // average number
double knr = 10.0; // spring constant of the number
double cvol[2] = {350, 350}; // average volume
double kvol = 10.0;
int nstvmov = 100;
int nequil = 100000;
int nsteps = 10000000;
double rho = 0.9;
double tp = 1.15;
double rcdef = 1e9; /* half-box cutoff */
double amp = 0.2; /* Monte Carlo move size */
const char *fnpos = "lj.pos";
int dopr = 0;



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

  return exp(-beta * de) * lj->vol / (n + 1);
}



/* volume move */
int lj_vmove(lj_t *lj, double beta, double cvol, double kvol)
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
     - kvol * (vn - vo) * (vn + vo - 2 * cvol);
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
    lj_energy(lj);
  }
  return acc;
}


/* particle move */
int lj_nmove(lj_t *lj, double beta, double cnr, double knr)
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
       + knr * (2 * cnr - 2 * n - 1);
    acc = ( dw > 0 || rand01() < exp(dw) );
    if ( acc ) {
      vcopy(lj->x[n], xi);
      // ignore v and f
      lj->n = n + 1;
      lj_setrho(lj, lj->n/lj->vol);
      lj->dof = lj->n * D;
      lj_energy(lj);
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
       - knr * (2 * cnr - 2 * n + 1);
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
      lj_energy(lj);
    }
  }

  return acc;
}


static void vexchange(lj_t *lj[2], double cvol[2],
    double beta, double time, av_t *avp)
{
  double p[2], dp, dvol, k;
  int i;

  for ( i = 0; i < 2; i++ )
    p[i] = av_getave(&avp[i]);
  dp = p[0] - p[1];
  k = ( lj[0]->n / lj[0]->vol / lj[0]->vol
      + lj[1]->n / lj[1]->vol / lj[1]->vol ) / beta;
  dvol = dp / k / time;

  /* adjust the volume exchange */
  if ( dvol > 0 ) {
    if ( dvol > cvol[1] * 0.5 )
      dvol = cvol[1] * 0.5;
  } else {
    if ( -dvol > cvol[0] * 0.5 )
      dvol = -cvol[0] * 0.5;
  }

  cvol[0] += dvol;
  cvol[1] -= dvol;
}


static void nexchange(lj_t *lj[2], double cnr[2],
    double time, av_t *avwid)
{
  double wid[2], dwid, dnr, k;
  int i;

  for ( i = 0; i < 2; i++ )
    wid[i] = av_getave(&avwid[i]);
  dwid = wid[0] - wid[1];
  k = -5.0 / (lj[0]->n + lj[1]->n);
  dnr = dwid / k / time;

  /* adjust the number exchange */
  if ( dnr > 0 ) {
    if ( dnr > cnr[1] * 0.5 )
      dnr = cnr[1] * 0.5;
  } else {
    if ( -dnr > cnr[0] * 0.5 )
      dnr = -cnr[0] * 0.5;
  }

  cnr[0] += dnr;
  cnr[1] -= dnr;
}


int main(void)
{
  int i, t, acc[2], nacc[2] = {0, 0}, vacc[2] = {0, 0};
  lj_t *lj[2];
  av_t avep[2], avacc[2], avnr[2], avvol[2];
  av_t avp[2], avwid[2], avgp[2], avgwid[2];
  double beta = 1/tp;

  mtscramble(time(NULL));

  for ( i = 0; i < 2; i++ ) {
    lj[i] = lj_open(nmax, rho, rcdef, dopr);
    lj_energy(lj[i]);
  }

  /* equilibration */
  for ( t = 1; t <= nequil; t++ ) {
    for ( i = 0; i < 2; i++ ) {
      lj_metro(lj[i], amp, beta);
    }
  }

  /* production */
  for ( i = 0; i < 2; i++ ) {
    av_clear(&avep[i]);
    av_clear(&avp[i]);
    av_clear(&avacc[i]);
    av_clear(&avnr[i]);
    av_clear(&avvol[i]);
    av_clear(&avwid[i]);
  }

  for ( t = 1; t <= nsteps; t++ ) {
    for ( i = 0; i < 2; i++ ) {
      acc[i] = lj_metro(lj[i], amp, beta);
      av_add(&avacc[i], acc[i]);
      av_add(&avp[i], lj_calcp(lj[i], tp));
      av_add(&avwid[i], lj_widom(lj[i], beta));
      if ( t % 1 == 0 ) {
        //nacc[i] = lj_nmove(lj[i], beta, cnr[i], knr);
        //av_add(&avnr[i], lj[i]->n);
        vacc[i] = lj_vmove(lj[i], beta, cvol[i], kvol);
        av_add(&avvol[i], lj[i]->vol);
      }
    }
    if ( t % nstvmov == 0 ) {
      vexchange(lj, cvol, beta, t / nstvmov, avp);
      //nexchange(lj, cnr, t / nstvmov, avwid);
      //printf("cvol %g, %g, vol %g, %g, p %g, %g, cnr %g, %g\n",
      //    cvol[0], cvol[1], lj[0]->vol, lj[1]->vol,
      //    av_getave(&avp[0]), av_getave(&avp[1]),
      //    cnr[0], cnr[1]);
      //getchar();
    }
    if ( t % 100000 == 0 ) {
      double mu[2], p[2], fr[2], vol[2];
      for ( i = 0; i < 2; i++ ) {
        mu[i] = -tp*log(av_getave(&avwid[i]));
        p[i] = av_getave(&avp[i]);
        vol[i] = av_getave(&avvol[i]);
        //printf("i %d, vol %g, %g\n", i, lj[i]->vol, av_getave(&avvol[i]));
        fr[i] = mu[i] - p[i]*vol[i]/lj[i]->n;
      }
      //getchar();

      printf("t %d, acc %g, %g, n %d, %d, nr %g, %g (%g, %g), mu %g, %g, nacc %d, %d,\n"
             "  vol %g, %g, av. vol %g, %g (%g, %g), rho %g, %g,"
             "  p %g, %g, fr %g, %g, vacc %d, %d\n",
         t, av_getave(&avacc[0]), av_getave(&avacc[1]),
         lj[0]->n, lj[1]->n, av_getave(&avnr[0]), av_getave(&avnr[1]),
         cnr[0], cnr[1], mu[0], mu[1], nacc[0], nacc[1],
         lj[0]->vol, lj[1]->vol, vol[0], vol[1], cvol[0], cvol[1],
         cnr[0]/vol[0], cnr[1]/vol[1],
         p[0], p[1], fr[0], fr[1], vacc[0], vacc[1]);
    }
    if ( t % nstvmov == 0 ) {
      av_clear(&avp[i]);
      //av_clear(&avwid[i]);
    }
  }
  for ( i = 0; i < 2; i++ ) {
    lj_close(lj[i]);
  }
  /*
  printf("tp %g, ep %g, nr %g, vol %g, acc %g%%, p %g, bmu %g\n",
      tp, av_getave(avep)/nr, av_getave(avp),
      av_getave(avnr), av_getave(avvol),
      100.*av_getave(avacc), av_getave(avp),
      -log(av_getave(avwid)*lj->vol/nr) );
  */
  return 0;
}

