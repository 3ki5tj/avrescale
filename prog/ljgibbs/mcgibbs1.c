/* pseudo-Gibbs ensemble Monte Carlo simulation */
#include "lj_npr.h"
#include "av.h"
#include "bpav.h"


int n = 108;
int nequil = 100000;
int nsteps = 50000000;
double rho[2] = {0.05, 0.6};
double tp = 1.15;
double rcdef = 1e9; /* half-box cutoff */
int nstvol = 100000;
int nstreport = 1000000;


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


static int volmove(lj_t *lj[2],
    double beta, double time,
    av_t *avbp, double gdbp[2],
    av_t *avw, double dlnw[2])
{
  double bp[2], vol[2], ln[2], s[2], ep6[2], ep12[2], eptail[2];
  double w[2], a[2][2], det;
  double delbp, dellnw, dvol[2], k, kref, de, ds, dw;
  int i, j, acc = 0;

  for ( i = 0; i < 2; i++ ) {
    bp[i] = av_getave(&avbp[i]);
    w[i] = av_getave(&avw[i]);
  } 
  delbp = bp[0] - bp[1];
  dellnw = log((w[0] + 1e-300)/(w[1] + 1e-300));

  /* solve the following equations
   * -delbp  = a00 * dvol0 + a01 * dvol1
   * -dellnw = a10 * dvol0 + a11 * dvol1
   */
  for ( i = 0; i < 2; i++ ) {
    if ( gdbp[i] > 0 ) {
      gdbp[i] = -lj[i]->n / lj[i]->vol / lj[i]->vol;
    }
  }

  a[0][0] = gdbp[0];
  a[0][1] = -gdbp[1];
  a[1][0] = dlnw[0];
  a[1][1] = -dlnw[1];
  det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
  dvol[0] = (-delbp * a[1][1] + dellnw * a[0][1]) / det;
  dvol[1] = ( delbp * a[1][0] - dellnw * a[0][0]) / det;

  /* modify the values by 1/t */
  for ( i = 0; i < 2; i++ )
    dvol[i] /= time;

  printf("vol %g%+g=%g, %g%+g=%g; rho %g, %g, (%g, %g, %g, %g) det %g, bp %g, %g, w %g, %g\n",
      lj[0]->vol, dvol[0], lj[0]->vol + dvol[0],
      lj[1]->vol, dvol[1], lj[1]->vol + dvol[1],
      lj[0]->n/lj[0]->vol, lj[1]->n/lj[1]->vol,
      a[0][0], a[0][1], a[1][0], a[1][1], det,
      bp[0], bp[1], w[0], w[1]);
  //getchar();

  vol[0] = lj[0]->vol + dvol[0];
  vol[1] = lj[1]->vol + dvol[1];
  if ( vol[0] > vol[1] && vol[1] > 0 ) { /* filter out unreasonable moves */
    /* compute the energy change */
    acc = 1;
    for ( i = 0; i < 2; i++ ) {
      ln[i] = pow(vol[i], 1./D);
      s[i] = vol[i] / lj[i]->vol;
      s[i] = 1 / (s[i] * s[i]);
      ep6[i] = lj[i]->ep6 * s[i];
      ep12[i] = lj[i]->ep12 * s[i] * s[i];
      eptail[i] = lj_gettail(ln[i]/2, lj[i]->n / vol[i], lj[i]->n, NULL);
      de = ((ep12[i] - ep6[i]) - (lj[i]->ep12 - lj[i]->ep6))
         + eptail[i] - lj[i]->epot_tail;
      ds = lj[i]->n * log( vol[i] / lj[i]->vol );
      dw = ds - beta * de;
      if ( dw < -10 ) {
        //printf("i %d, dw %g, ds %g, de %g\n", i, dw, ds, de);
        acc = 0;
        break;
      }
    }

    if ( acc ) {
      for ( i = 0; i < 2; i++ )
        lj_changevol(lj[i], vol[i], ep6[i], ep12[i]);
    }
  }

/*
  printf("t %g, vol %g, %g, dvol %g, %g, acc %d, ds %g, de %g, dw %g, "
         "bp %g, %g, delbp %g, k %g (k*n %g)\n",
         time, lj[0]->vol, lj[1]->vol, dvol[0], dvol[1], acc,
         ds, de, dw, bp[0], bp[1], delbp, k, k*lj[0]->n);
*/

  if ( acc ) {
    for ( i = 0; i < 2; i++ ) {
      av_clear(&avbp[i]);
      av_clear(&avw[i]);
    }
  }
  return acc;
}


int main(void)
{
  int i, t, acc[2], volacc;
  lj_t *lj[2];
  av_t avbp[2], avw[2], avacc[2], avvol[2], avvolacc[1];
  bpav_t bpav[2];
  double bp[2], beta = 1.0/tp, amp, gbp[2], gdbp[2], w[2], dlnw[2];

  mtscramble(time(NULL));
  for ( i = 0; i < 2; i++ ) {
    lj[i] = lj_open(n, rho[i], rcdef);
    lj[i]->dof = n * D;
    lj_energy(lj[i]);
  }

  /* equilibration */
  for ( t = 1; t <= nequil; t++ ) {
    for ( i = 0; i < 2; i++ ) {
      lj_metro(lj[i], 0.1/rho[i], 1/tp);
    }
  }

  /* production */
  for ( i = 0; i < 2; i++ ) {
    av_clear(&avbp[i]);
    av_clear(&avw[i]);
    av_clear(&avacc[i]);
    bpav_clear(&bpav[i]);
    av_clear(&avvol[i]);
  }
  av_clear(avvolacc);

  for ( t = 1; t <= nsteps; t++ ) {
    for ( i = 0; i < 2; i++ ) {
      double rho = lj[i]->n / lj[i]->vol;
      // empirical formula for the MC move size
      amp = 0.1 / rho;
      acc[i] = lj_metro(lj[i], amp, beta);
      av_add(&avacc[i], acc[i]);
      bpav_add(&bpav[i], lj[i], beta);
      av_add(&avbp[i], bpav[i].bp);
      av_add(&avw[i], bpav[i].w * lj[i]->vol / (lj[i]->n + 1));
      //printf("i %d, bp %g, %g\n", i, bp[i], av_getave(&avgp[i])); getchar();
    }

    if ( t % nstvol == 0 ) {
      for ( i = 0; i < 2; i++ ) {
        gbp[i] = bpav_get(&bpav[i], &gdbp[i], &w[i], &dlnw[i]);
        w[i] *= lj[i]->vol / (lj[i]->n + 1);
        dlnw[i] += 1./lj[i]->vol;
      }
      volacc = volmove(lj, beta, t / nstvol, avbp, gdbp, avw, dlnw);
      av_add(avvolacc, volacc);
      for ( i = 0; i < 2; i++ ) {
        av_add(&avvol[i], lj[i]->vol);
      }
    }


    if ( t % nstreport == 0 ) {
      for ( i = 0; i < 2; i++ ) {
        gbp[i] = bpav_get(&bpav[i], &gdbp[i], &w[i], &dlnw[i]);
        w[i] *= lj[i]->vol / (lj[i]->n + 1);
        dlnw[i] += 1./lj[i]->vol;
        //bpav_clear(&bpav[i]);
      }
      printf("t %d, rho %g %g, bp %g %g, dbp %g, %g, w %g, %g, dlnw %g, %g, "
             "vmove %g%%, acc %g%%, %g%%\n",
          t, lj[0]->n/lj[0]->vol, lj[1]->n/lj[1]->vol,
          gbp[0], gbp[1], gdbp[0], gdbp[1],
          w[0], w[1], dlnw[0], dlnw[1],
          av_getave(avvolacc)*100,
          av_getave(&avacc[0])*100, av_getave(&avacc[1])*100);
      //getchar();
    }
  }

  for ( i = 0; i < 2; i++ ) {
    lj_close(lj[i]);
  }
  return 0;
}

