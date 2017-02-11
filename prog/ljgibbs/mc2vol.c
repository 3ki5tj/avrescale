/* pseudo-Gibbs ensemble Monte Carlo simulation */
#include <time.h>
#include "lj.h"
#include "av.h"
#include "bpav.h"


int n[2] = {200, 200};
int nequil = 1000000;
int nsteps = 100000000;
double den[2] = {0.13, 0.43};
double rholimit[2] = {0.25, 0.35};
double tp = 1.25;
double rcdef = 1e9; /* half-box cutoff */
int dopr = 0;
int nstvol = 1000;
int nstreport = 100000;
double mageql = 0.01; /* scaling magnitude during equilibration */
int flipmin = 100; /* minimal number of flips to start stage 1 */
int force_simple = 0; /* always use ideal-gas formula */


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
    av_t *avbp, double gdbp[2], av_t *avwv, double gwv[2], double gdlnwv[2])
{
  double bp[2], nvol[2], ln[2], s[2], ep6[2], ep12[2], eptail[2];
  double rho[2], wv[2];
  double a[2][2] = {{0, 0}, {0, 0}}, det = 0, mag;
  double delbp, dellnwv, dlnvol[2], de, ds, dse;
  int i, acc = 0, holdv1;

  for ( i = 0; i < 2; i++ ) {
    bp[i] = av_getave(&avbp[i]);
    wv[i] = av_getave(&avwv[i]);
  }
  delbp = bp[0] - bp[1];
  /* approximate formula for ln(wv[0] / wv[1]) */
  dellnwv = (wv[0] - wv[1]) * 2 / (gwv[0] + gwv[1]);
  /* update the number of flips */
  if ( fabs(del[0]) > 0 && fabs(del[1]) ) {
    if ( delbp * del[0] < 0 ) nflips[0] += 1;
    if ( dellnwv * del[1] < 0 ) nflips[1] += 1;
  }
  del[0] = delbp;
  del[1] = dellnwv; /* difference of < exp(-dU) V / (N+1) > */

  /* compute the current density and specific volume */
  for ( i = 0; i < 2; i++ ) {
    rho[i] = lj[i]->n / lj[i]->vol;
  }

  holdv1 = (gdbp[1] < 0);

  if ( simple ) {
    /* Use ideal-gas formula
     * pressure depends more critically on the liquid-phase
     * -delbp  = (dbp0/dlnV0) dlnV0 - (dbp1/dlnV1) dlnV1
     *         ~ -rho0 dlnV0 + rho1 dlnV1
     * -dellnwv  = (dlnw0/dlnV0) dlnV0 - (dlnw1/dlnV1) dlnV1
     *           ~ dlnV0 - dlnV1 */
    a[0][0] = -rho[0];
    a[0][1] =  rho[1];
    a[1][0] =  1;
    a[1][1] = -1;
  } else {
    /* safeguarding the derivatives */
    for ( i = 0; i < 2; i++ ) {
      if ( gdbp[i] < 0.1 ) gdbp[i] = 0.1;
      if ( gdlnwv[i] < 0.1 ) gdlnwv[i] = 0.1;
    }

    /* solve the following equations
     * -delbp   = a00 * dlnvol0 + a01 * dlnvol1
     * -dellnwv = a10 * dlnvol0 + a11 * dlnvol1
     */
    a[0][0] = -gdbp[0] * rho[0]; /*  dbp0/dlnvol0 = -dbp0/drho0*rho0 */
    a[0][1] =  gdbp[1] * rho[1]; /* -dbp1/dlnvol1 =  dbp1/drho1*rho1 */
    /* by thermodynamics gdlnwv == gdbp
     * i.e. rho * d(mu)/d(rho) = d p/d(rho) */
    // a[1][0] =  gdlnwv[0]; /*  dlnwv0/dlnvol0 */
    // a[1][1] = -gdlnwv[1]; /* -dlnwv1/dlnvol1 */
    a[1][0] =  gdbp[0]; /*  dlnwv0/dlnvol0 =  dbp0/drho0 */
    a[1][1] = -gdbp[1]; /* -dlnwv1/dlnvol1 = -dbp1/drho1 */
  }
  det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
  dlnvol[0] = (-delbp * a[1][1] + dellnwv * a[0][1]) / det;
  dlnvol[1] = ( delbp * a[1][0] - dellnwv * a[0][0]) / det;
  /* volume of the liquid phase too small */
  if ( holdv1 && dlnvol[1] < 0.1 ) dlnvol[1] = 0.1;
  //printf("%g, %g, %g, %g; %g %g; %g, %g; mag %g\n", a[0][0], a[0][1], a[1][0], a[1][1], delbp, dellnwv, dlnvol[0], dlnvol[1], mag0);

  /* modify the values by 1/t */
  if ( ntimes != NULL ) {
    mag = 1. / (*ntimes + 2.0/mag0);
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
           "(%.4f, %8.3f, %8.3f, %8.3f) det %8.2f, dbp %+.4f%+.4f=%+8.3f, dwv %6.3f%+6.3f=%+8.3f, mag %g\n",
        lj[0]->vol, exp(dlnvol[0]), nvol[0],
        lj[1]->vol, exp(dlnvol[1]), nvol[1],
        rho[0], rho[1],
        a[0][0], a[0][1], a[1][0], a[1][1], det,
        bp[0], bp[1], del[0], wv[0], wv[1], del[1], mag);
    //getchar();
  }

  if ( lj[0]->n/nvol[0] < rholimit[0]
    && lj[1]->n/nvol[1] > rholimit[1] ) { /* filter out unreasonable moves */
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
  av_t avbp[2], avwv[2], avacc[2], avvolacc[1];
  bpav_t bpav[2];
  double beta = 1.0/tp, amp, mag;
  double rho[2], spv[2];
  double gbp[2], gdbp[2], gw[2], gdlnw[2], gwv[2], gdlnwv[2], del[2] = {0, 0};

  //mtscramble(time(NULL));
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
  }
  av_clear(avvolacc);

  for ( t = 1; stage < 2 || t <= nsteps; t++ ) {
    for ( i = 0; i < 2; i++ ) {
      rho[i] = lj[i]->n / lj[i]->vol;
      spv[i] = lj[i]->vol / (lj[i]->n + 1);
      /* recalibrate the energy to avoid accumulated
       * numerical error, unnecessary but only for safety */
      if ( t % 10000 == 0 ) lj_energy(lj[i]);
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
        gbp[i] = bpav_get(&bpav[i], lj[i], &gdbp[i], &gw[i], &gdlnw[i]);
        gwv[i] = gw[i] * spv[i];
        gdlnwv[i] = 1 + gdlnw[i];
      }

      if ( stage == 0 ) { /* equilibration */
        mag = mageql;
        pvol = NULL;
        simple = 1; /* use simplified formula in the equilibration phase */
        verbose = 1;
      } else if ( stage == 1 ) {
        mag = mageql;
        pvol = &nvol;
        simple = 1;
#if 0
        mag = 0; /* stop updating */
        pvol = NULL;
        simple = 0;
#endif
        verbose = 1;
      } else { /* production */
        mag = mageql;
        pvol = &nvol;
        simple = force_simple;
        verbose = (t > 0 && t % nstreport == 0);
      }

      /* adjust the two volumes */
      volacc = volmove(lj, beta, mag, pvol,
          simple, del, nflips, verbose, avbp, gdbp, avwv, gwv, gdlnwv);

      /* check if we need to change stages */
      if ( ( stage == 0 && nflips[0] >= flipmin && nflips[1] >= flipmin ) 
        || ( stage == 1 && t >= nequil && gdbp[0] > 0 && gdbp[1] > 0 ) ) {
        stage += 1;
        fprintf(stderr, "starting stage %d at t %d, rho %g, %g, nflips %d, %d, gdbp %g, %g\n",
            stage, t, rho[0], rho[1], nflips[0], nflips[1], gdbp[0], gdbp[1]);
        t = 0;
        nflips[0] = nflips[1] = 0;
        if ( stage == 1 ) {
          for ( i = 0; i < 2; i++ ) bpav_clear(&bpav[i]);
        }
      }

      if ( stage == 2 ) {
        /* production stage */
        av_add(avvolacc, volacc);
      }
    }


    if ( t > 0 && t % nstreport == 0 ) {
      /* get the global averages */
      for ( i = 0; i < 2; i++ ) {
        gbp[i] = bpav_get(&bpav[i], lj[i], &gdbp[i], &gw[i], &gdlnw[i]);
        gwv[i] = gw[i] * spv[i];
        gdlnwv[i] = 1 + gdlnw[i];
        //bpav_clear(&bpav[i]);
      }
      printf("t %d@%d: rho %8.6f %8.6f, bp %8.6f %8.6f, dbp %g %g, wv %7.4f %7.4f, dlnwv %g %g, "
             "vmove %.3f%%, acc %.0f%% %.0f%%, flips %d %d, del %+.3f %+6.2f\n",
          stage, t, rho[0], rho[1],
          gbp[0], gbp[1], gdbp[0], gdbp[1],
          gwv[0], gwv[1], gdlnwv[0], gdlnwv[1],
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

