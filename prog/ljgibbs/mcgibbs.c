/* pseudo-Gibbs ensemble Monte Carlo simulation */
#include "lj.h"
#include "av.h"


int n = 108;
int nequil = 100000;
int nsteps = 50000000;
double rho = 0.3;
double tp = 1.15;
double rcdef = 1e9; /* half-box cutoff */
double amp = 0.2; /* Monte Carlo move size */
const char *fnpos1 = "lj0.pos";
const char *fnpos2 = "lj1.pos";
double lnvamp = 0.05;
double pres = 1.55; /* corresponding to rho */
int nstvmov = 1000;
int nstreport = 10000;


static void lj_vscale(lj_t *lj, double vol)
{
  double s = pow(vol/lj->vol, 1/D);
  int i;
  lj->vol = vol;
  lj_setrho(lj, lj->n / vol);
  for ( i = 0; i < lj->n; i++ )
    vsmul(lj->x[i], s);
  lj_force(lj);
}

static void vexchange(lj_t *lj[2],
    double tp, double time, av_t *avp)
{
  double p0, p1, dp, dvol, vol0, vol1, k;
  p0 = av_getave(&avp[0]);
  p1 = av_getave(&avp[1]);
  dp = p0 - p1;
  k = (lj[0]->n*tp/lj[0]->vol/lj[0]->vol + lj[1]->n*tp/lj[1]->vol/lj[1]->vol);
  /* if system 1 has a higher pressure than system 2
   * expand system 1, otherwise expand system 2 */
  dvol = dp / k / time;
  vol0 = lj[0]->vol + dvol;
  if ( vol0 < 0 ) {
    vol0 = lj[0]->vol * 0.5;
    dvol = vol0 - lj[0]->vol;
  }
  vol1 = lj[1]->vol - dvol;
  if ( vol1 < 0 ) {
    vol1 = lj[1]->vol * 0.5;
    dvol = lj[1]->vol - vol1;
    vol0 = lj[0]->vol + dvol;
  }
  lj_vscale(lj[0], vol0);
  lj_vscale(lj[1], vol1);
  printf("t %g, p %g, %g, dp %g, k %g (k*n %g), vol %g, %g, dvol %g\n",
      time, p0, p1, dp, k, k*lj[0]->n, lj[0]->vol, lj[1]->vol, dvol);
  av_clear(&avp[0]);
  av_clear(&avp[1]);
}

int main(void)
{
  int i, t, acc[2];
  lj_t *lj[2];
  av_t avp[2], avacc[2], avgp[2], avvol[2];
  double p[2];

  mtscramble(time(NULL));
  for ( i = 0; i < 2; i++ ) {
    lj[i] = lj_open(n, rho, rcdef);
    lj_energy(lj[i]);
  }

  /* equilibration */
  for ( t = 1; t <= nequil; t++ ) {
    for ( i = 0; i < 2; i++ ) {
      lj_metro(lj[i], amp, 1/tp);
    }
  }

  /* production */
  for ( i = 0; i < 2; i++ ) {
    av_clear(&avp[i]);
    av_clear(&avacc[i]);
    av_clear(&avgp[i]);
    av_clear(&avvol[i]);
  }

  for ( t = 1; t <= nsteps; t++ ) {
    for ( i = 0; i < 2; i++ ) {
      acc[i] = lj_metro(lj[i], amp, 1/tp);
      av_add(&avacc[i], acc[i]);
      p[i] = lj_calcp(lj[i], tp);
      av_add(&avp[i], p[i]);
      av_add(&avgp[i], p[i]);
      //printf("i %d, p %g, %g\n", i, p[i], av_getave(&avgp[i])); getchar();
    }

    if ( t % nstvmov == 0 ) {
      vexchange(lj, tp, t / nstvmov, avp);
      for ( i = 0; i < 2; i++ ) {
        av_add(&avvol[i], lj[i]->vol);
      }
    }

    if ( t % nstreport == 0 ) {
      printf("t %d, rho %g %g, p %g %g\n",
          t, lj[0]->n/av_getave(&avvol[0]), lj[1]->n/av_getave(&avvol[1]),
          av_getave(&avgp[0]), av_getave(&avgp[1]));
    }
  }
  for ( i = 0; i < 2; i++ ) {
    lj_close(lj[i]);
  }
  return 0;
}

