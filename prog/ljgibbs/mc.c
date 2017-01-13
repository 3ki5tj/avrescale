/* basic Monte Carlo simulation */
#include "lj_npr.h"
#include "av.h"
#include "bpav.h"



int n = 108;
int nequil = 100000;
int nsteps = 5000000;
double rho = 0.605; // 0.07; // 0.605;
//double rho = 0.07;
double tp = 1.15;
double rcdef = 1e9; /* half-box cutoff */
const char *fnpos = "lj.pos";



int main(void)
{
  int t, acc;
  lj_t *lj;
  double epsm = 0, accsm = 0;
  double beta = 1/tp, bp, dbp, w, dlnw, amp = 0.1/rho;
  bpav_t bpav[1];

  mtscramble(time(NULL));
  bpav_clear(bpav);
  lj = lj_open(n, rho, rcdef);
  lj->dof = n * D;
  //lj->vol -= 0.5; lj_setrho(lj, n / lj->vol); 
  //lj->vol += 0.5; lj_setrho(lj, n / lj->vol); 
  lj_energy(lj);
  /* equilibration */
  for ( t = 1; t <= nequil; t++ )
    lj_metro(lj, amp, 1/tp);
  /* production */
  for ( t = 1; t <= nsteps; t++ ) {
    acc = lj_metro(lj, amp, 1/tp);
    epsm += lj->epot;
    accsm += acc;
    bpav_add(bpav, lj, beta);
  }
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
  bp = bpav_get(bpav, &dbp, &w, &dlnw);
  w *= lj->vol / lj->n;
  dlnw += 1./lj->vol;
  printf("rho %g, tp %g, vol %g, ep %g, bp %g, dbp %g, w %g, lnw %g, dlnw %g, acc %g%%\n",
      rho, tp, lj->vol, epsm / nsteps, bp, dbp, w, log(w), dlnw,
      100.*accsm / nsteps);
  return 0;
}

