/* solving equal pressure lines */
#include "ljeos.h"

#define N 10000
double rho[N], p[N], mu[N], rhop[N], rhomu[N];

int main(void)
{
  int i, j, imin = 0;
  double tp = 1.15, f, rhog = 0, rhol = 0;
  for ( i = 0; i < N; i++ ) {
    rho[i] = 1.0 * (i + 1) / N;
    ljeos3d_get(rho[i], tp, &p[i], &f, &mu[i]);
    mu[i] += tp * log(rho[i]);
    //printf("rho %g, p %g, mu %g\n", rho[i], p[i], mu[i]);
  }
  //exit(1);
  for ( i = 0; i < N; i++ ) {
    // for each density, find another density
    // with the same pressure
    for ( j = N - 1; j >= i; j-- ) {
      if ( p[j] < p[i] ) {
        rhop[i] = rho[j];
        break;
      }
    }
    for ( j = N - 1; j >= i; j-- ) {
      if ( mu[j] < mu[i] ) {
        rhomu[i] = rho[j];
        break;
      }
    }
    if ( j >= i && imin == 0 ) {
      imin = i;
    }
    // stop searching once p starts decreasing
    if ( p[i+1] <= p[i] ) break;
  }

  for ( j = imin; j < i; j++ ) {
    if ( fabs(rhop[j] - rhomu[j]) < 2./N ) {
      rhog = rho[j];
      rhol = rhop[j];
    }
    printf("%g %g %g %g %g\n", rho[j], rhop[j], rhomu[j], p[j], mu[j]);
  }
  printf("# gas %g, liquid %g\n", rhog, rhol);
  return 0;
}

