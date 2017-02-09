/* pseudo-Gibbs ensemble Monte Carlo simulation */
#include <time.h>
#include "lj.h"
#include "av.h"
#include "bpav.h"


int n[2] = {100, 100};
int nequil = 500000;
int nsteps = 50000000;
double den[2] = {0.05, 0.6};
double tp = 1.15;
double rcdef = 1e9; /* half-box cutoff */
int dopr = 0;
int nstvol = 1000;
int nstreport = 1000000;
double mageql = 0.1; /* scaling magnitude during equilibration */
int flipmin = 10; /* minimal number of flips to start stage 1 */


/* change the volume */
static void lj_changevol(lj_t *lj, double vol, double ep6, double ep12)
{
  int j, n = lj->n;
  double s = pow(vol / lj->vol, 1.0/D);
  lj_setrho(lj, n / vol);
  for ( j = 0; j < n; j++ )
    vsmul(lj->x[j], s);
  // to check the energy change rigorously
  //lj_energy(lj);
  //printf("i %d, ep6 %g, %g, ep12 %g, %g, ep0 %g, %g, ep %g, %g, vir %g, %g\n",
  //    i, lj->ep6, ep6, lj->ep12, ep12, lj->ep0, ep12-ep6,
  //    lj->epot, ep12-ep6+eptail, lj->vir, 12*ep12-6*ep6);
  //getchar();
  /* update the energies by scaling
   * not so efficient if we have to update r2ij */
  lj->ep6 = ep6;
  lj->ep12 = ep12;
  lj->ep0 = ep12 - ep6;
  lj->epot = lj->ep0 + lj->epot_tail;
  lj->vir = 12 * ep12 - 6 * ep6;
  //lj_energy(lj);
}


/* volume move
 * adjust the volumes of the boxes to make their pressure
 * and chemical potential match
 * `simple`: use simplified formula for derivatives
 * if magnitude `mag` is equal 0, then we will use the 1/t scheme */
static int volmove(lj_t *lj[2], double beta, double mag0, int *ntimes,
    int simple, double del[2], int nflips[2], int verbose,
    av_t *avbp, double gdbp[2], av_t *avwv, double gdwv[2])
{
  double bp[2], nvol[2], ln[2], s[2], ep6[2], ep12[2], eptail[2];
  double rho[2], spv[2], wv[2];
  double a[2][2] = {{0, 0}, {0, 0}}, det = 0, mag;
  double delbp, delwv, dlnvol[2], de, ds, dse;
  int i, acc = 0;

  for ( i = 0; i < 2; i++ ) {
    bp[i] = av_getave(&avbp[i]);
    wv[i] = av_getave(&avwv[i]);
  }
  delbp = bp[0] - bp[1];
  delwv = wv[0] - wv[1];
  /* update the number of flips */
  if ( fabs(del[0]) > 0 && fabs(del[1]) ) {
    if ( delbp * del[0] < 0 ) nflips[0] += 1;
    if ( delwv * del[1] < 0 ) nflips[1] += 1;
  }
  del[0] = delbp;
  del[1] = delwv; /* difference of < exp(-dU) V / (N+1) > */

  for ( i = 0; i < 2; i++ ) {
    rho[i] = lj[i]->n / lj[i]->vol;
    spv[i] = lj[i]->vol / (lj[i]->n + 1);
  }

  if ( simple ) {
    /* Use ideal-gas formula
     * pressure depends more critically on the liquid-phase
     * -delbp  = (dbp0/dlnV0) dlnV0 - (dbp1/dlnV1) dlnV1
     *         ~ -rho0 dlnV0 + rho1 dlnV1
     * -delwv  = (dw0/dlnV0) dlnV0 - (dw1/dlnV1) dlnV1
     *         ~ (1/rho0) dlnV0 - (1/rho1) dlnV1 */
    a[0][0] = -rho[0];
    a[0][1] =  rho[1];
    a[1][0] =  spv[0];
    a[1][1] = -spv[1];
  } else {
    /* use ideal gas values for unstable values */
    for ( i = 0; i < 2; i++ ) {
      if ( gdbp[i] < 0.1 ) gdbp[i] = 0.1;
      if ( gdwv[i] < 0.1*spv[i] ) gdwv[i] = 0.1*spv[i];
    }

    /* solve the following equations
     * -delbp  = a00 * dlnvol0 + a01 * dlnvol1
     * -delwv  = a10 * dlnvol0 + a11 * dlnvol1
     */
    a[0][0] = -gdbp[0] * rho[0]; /*  dbp0/dlnvol0 = -dbp0/drho0*rho0 */
    a[0][1] =  gdbp[1] * rho[1]; /* -dbp1/dlnvol1 =  dbp1/drho1*rho1 */
    a[1][0] =  gdwv[0]; /*  dwv0/dlnvol0 */
    a[1][1] = -gdwv[1]; /* -dwv1/dlnvol1 */
  }
  det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
  dlnvol[0] = (-delbp * a[1][1] + delwv * a[0][1]) / det;
  dlnvol[1] = ( delbp * a[1][0] - delwv * a[0][0]) / det;
  //printf("%g, %g, %g, %g; %g %g; %g, %g; mag %g\n", a[0][0], a[0][1], a[1][0], a[1][1], delbp, delw, dlnvol[0], dlnvol[1], mag0);

  /* modify the values by 1/t */
  if ( ntimes != NULL ) {
    mag = 1. / (*ntimes + 1);
    *ntimes += 1;
  } else {
    mag = mag0;
  }

  /* compute the new volumes */
  for ( i = 0; i < 2; i++ ) {
    dlnvol[i] *= mag;
    nvol[i] = lj[i]->vol * exp(dlnvol[i]);
  }

  if ( verbose ) {
    printf("vol %7.2f*%5.3f=%7.2f, %6.2f*%5.3f=%6.2f; rho %.4f, %.4f, "
           "(%.4f, %8.3f, %8.3f, %8.3f) det %8.2f, dbp %+.4f%+.4f=%+8.3f, dwv %6.3f%+6.3f=%+8.3f\n",
        lj[0]->vol, exp(dlnvol[0]), nvol[0],
        lj[1]->vol, exp(dlnvol[1]), nvol[1],
        rho[0], rho[1],
        a[0][0], a[0][1], a[1][0], a[1][1], det,
        bp[0], bp[1], del[0], wv[0], wv[1], del[1]);
    //getchar();
  }

  if ( lj[0]->n/nvol[0] < lj[1]->n/nvol[1] ) { /* filter out unreasonable moves */
    /* compute the energy change */
    acc = 1;
    for ( i = 0; i < 2; i++ ) {
      ln[i] = pow(nvol[i], 1./D);
      s[i] = nvol[i] / lj[i]->vol;
      s[i] = 1 / (s[i] * s[i]);
      ep6[i] = lj[i]->ep6 * s[i];
      ep12[i] = lj[i]->ep12 * s[i] * s[i];
      eptail[i] = lj_gettail(ln[i]/2, lj[i]->n / nvol[i], lj[i]->n, NULL);
      de = ((ep12[i] - ep6[i]) - (lj[i]->ep12 - lj[i]->ep6))
         + eptail[i] - lj[i]->epot_tail;
      ds = lj[i]->n * log( nvol[i] / lj[i]->vol );
      dse = ds - beta * de;
      if ( dse < -10 ) {
        //printf("i %d, dse %g, ds %g, de %g\n", i, dse, ds, de);
        acc = 0;
        break;
      }
    }

    if ( acc ) {
      /* commit the volume changes */
      for ( i = 0; i < 2; i++ )
        lj_changevol(lj[i], nvol[i], ep6[i], ep12[i]);
    }
  }

  if ( acc ) {
    /* clear the block averages
     * only need to do so if the move is accepted */
    for ( i = 0; i < 2; i++ ) {
      av_clear(&avbp[i]);
      av_clear(&avwv[i]);
    }
  }
  return acc;
}


int main(void)
{
  int i, t, acc[2], volacc, nvol = 0, nflips[2] = {0, 0};
  int stage = 0, simple, *pvol = NULL, verbose;
  lj_t *lj[2];
  av_t avbp[2], avwv[2], avacc[2], avvol[2], avvolacc[1];
  bpav_t bpav[2];
  double beta = 1.0/tp, amp, mag;
  double rho[2], spv[2];
  double gbp[2], gdbp[2], gw[2], gdw[2], gdwv[2], del[2] = {0, 0};

  mtscramble(time(NULL));
  for ( i = 0; i < 2; i++ ) {
    lj[i] = lj_open(n[i], den[i], rcdef, dopr);
    lj[i]->dof = n[i] * D; /* MC */
    lj_energy(lj[i]);
  }

  /* equilibration */
  for ( t = 1; t <= nequil; t++ ) {
    for ( i = 0; i < 2; i++ ) {
      amp = 0.1 / den[i]; /* empirical MC move size */
      lj_metro(lj[i], amp, 1/tp);
    }
  }

  /* production */
  for ( i = 0; i < 2; i++ ) {
    av_clear(&avbp[i]);
    av_clear(&avwv[i]);
    av_clear(&avacc[i]);
    bpav_clear(&bpav[i]);
    av_clear(&avvol[i]);
  }
  av_clear(avvolacc);

  for ( t = 1; stage < 2 || t <= nsteps; t++ ) {
    for ( i = 0; i < 2; i++ ) {
      rho[i] = lj[i]->n / lj[i]->vol;
      spv[i] = lj[i]->vol / (lj[i]->n + 1);
      amp = 0.1 / rho[i]; /* empirical MC move size */
      acc[i] = lj_metro(lj[i], amp, beta);
      av_add(&avacc[i], acc[i]);
      bpav_add(&bpav[i], lj[i], beta); /* global averages */
      /* two block averagers for the discrepencies in
       * beta * p and exp(-dU) * V / (N + 1) of the two boxes */
      av_add(&avbp[i], bpav[i].bp);
      av_add(&avwv[i], bpav[i].w * spv[i]);
      //if ( stage == 1) { printf("i %d, bp %g, %g\n", i, bpav[i].bp, av_getave(&avbp[i])); getchar(); }
    }

    if ( t % nstvol == 0 ) {
      /* get the global averages */
      for ( i = 0; i < 2; i++ ) {
        gbp[i] = bpav_get(&bpav[i], lj[i], &gdbp[i], &gw[i], &gdw[i]);
        gdwv[i] = spv[i] * (gw[i] + gdw[i]);
      }

      if ( stage == 0 ) { /* equilibration */
        mag = mageql;
        pvol = NULL;
        simple = 1; /* use simplified formula in the equilibration phase */
        verbose = 1;
      } else if ( stage == 1 ) {
        mag = 0; /* stop updating */
        pvol = NULL;
        simple = 0;
        verbose = 1;
      } else { /* production */
        mag = mageql;
        pvol = &nvol;
        simple = 0;
        verbose = (t > 0 && t % nstreport == 0);
      }

      /* adjust the two volumes */
      volacc = volmove(lj, beta, mag, pvol,
          simple, del, nflips, verbose, avbp, gdbp, avwv, gdwv);

      if ( stage == 0 && nflips[0] >= flipmin && nflips[1] >= flipmin ) {
        stage = 1;
        fprintf(stderr, "starting stage 1 at t %d, rho %g, %g, nflips %d, %d\n",
            t, rho[0], rho[1], nflips[0], nflips[1]);
        t = 0;
        nflips[0] = nflips[1] = 0;
        for ( i = 0; i < 2; i++ ) bpav_clear(&bpav[i]);
      } else if ( stage == 1 && t >= nequil ) {
        stage = 2;
        fprintf(stderr, "starting stage 2 at t %d, rho %g, %g, nflips %d, %d\n",
            t, rho[0], rho[1], nflips[0], nflips[1]);
        t = 0;
        nflips[0] = nflips[1] = 0;
      } else if ( stage == 2 ) {
        /* production stage */
        av_add(avvolacc, volacc);
        for ( i = 0; i < 2; i++ ) {
          av_add(&avvol[i], lj[i]->vol);
        }
      }
    }


    if ( t > 0 && t % nstreport == 0 ) {
      /* get the global averages */
      for ( i = 0; i < 2; i++ ) {
        gbp[i] = bpav_get(&bpav[i], lj[i], &gdbp[i], &gw[i], &gdw[i]);
        gdwv[i] = spv[i] * (gw[i] + gdw[i]);
        //bpav_clear(&bpav[i]);
      }
      printf("t %d@%d: rho %8.6f %8.6f, bp %8.6f %8.6f, dbp %g, %g, wv %7.4f, %7.4f, dwv %g, %g, "
             "vmove %.2f%%, acc %.0f%%, %.0f%%, flips %d, %d, del %+.3f, %+6.2f\n",
          stage, t, rho[0], rho[1],
          gbp[0], gbp[1], gdbp[0], gdbp[1],
          gw[0] * spv[0], gw[1] * spv[1], gdwv[0], gdwv[1],
          av_getave(avvolacc)*100,
          av_getave(&avacc[0])*100, av_getave(&avacc[1])*100,
          nflips[0], nflips[1], del[0], del[1]);
      //getchar();
    }
  }

  for ( i = 0; i < 2; i++ ) {
    lj_close(lj[i]);
  }
  return 0;
}

