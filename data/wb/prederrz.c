/* compute the analytical error */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int dof = 4785;
double kB = 0.001987191;
double tp = 300;
double gam = 0.33661; // from adp/avrescale.dat
int surf = 0; // surface
const char *fncorr = "fix/ene0fix_k.acf";
//const char *fncorr = "fix/ene0fix_invk.acf";
double initstd = 100;
double tottime = 1000000;

/* load the autocorrelation function from file */
static double *loadcorr(const char *fn, double *dt, int *cnt)
{
  double c;
  FILE *fp;
  char buf[1024];
  int i, ifr, n = 0, m = 0;
  double *xx, y, ave = 0, var = 0;

  if ( surf ) {
    c = kB*tp*kB*tp*(dof*0.5-1)*(dof*0.5+gam-2)/gam;
  }
  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot read %s\n", fn);
    return NULL;
  }
  fgets(buf, sizeof buf, fp);
  i = sscanf(buf + 1, "%d%lf%lf%lf%d", &n, &ave, &var, dt, &m);
  if ( i != 5 ) {
    fprintf(stderr, "cannot read 5 variables from %s\n", buf);
    return NULL;
  }
  var *= c * c;
  *cnt = m;
  xx = calloc(m + 1, sizeof(*xx));
  for ( i = 0; i <= m; i++ ) {
    fgets(buf, sizeof buf, fp);
    sscanf(buf, "%d%lf", &ifr, &y);
    xx[i] = y * var;
  }
  fprintf(stderr, "ave %g, var %g, dt %g, m %d\n", ave, var, *dt, m);
  fclose(fp);
  return xx;
}

static double geterrz(double z, const double *xx, double cnt,
    double dt, double time)
{
  int i, j, j0, tmax = (int) (time / dt + 0.5);
  double ss = 0, mul;
  double *du = calloc(tmax + 1, sizeof(double));

  for ( i = 1; i <= tmax; i++ )
    du[i] = z/(i*dt)*pow(i*dt/time, z);

  for ( i = 1; i <= tmax; i++ ) {
    for ( j = (i > cnt) ? i - cnt : 1; j <= i; j++ ) {
      mul = ( j == i ) ? 0.5 : 1.0;
      ss += du[i] * du[j] * xx[i-j] * mul;
    }
    //printf("i %d, udi %g, ss %g\n", i, udi, ss); getchar();
  }
  ss *= 2 * dt*dt;

  free(du);
  return ss;
}


/* white noise approximation */
static double geterrz_white(double z, const double *xx, double cnt,
    double dt, double time)
{
  int i, tmax = (int) (time / dt + 0.5);
  double ss, gamma = xx[0];
  for ( i = 1; i <= cnt; i++ ) gamma += 2 * xx[i];
  gamma *= dt;
  if ( fabs(2*z - 1) < 1e-10 ) {
    ss = gamma * z*z*log(time)/time;
  } else {
    ss = gamma * z*z/(2*z-1) * (1/time - 1/pow(time, 2*z));
  }
  return ss;
}


int main(int argc, char **argv)
{
  int cnt;
  double *xx, dt, z, zstep = 0.1, err, err_white, err_res;

  if ( argc > 1 ) tottime = atof( argv[1] );
  if ( argc > 2 ) zstep = atof( argv[2] );
  xx = loadcorr(fncorr, &dt, &cnt);
  if ( xx == NULL ) return -1;
  for ( z = 0.1; z <= 10; z += zstep ) {
    err = geterrz(z, xx, cnt, dt, tottime);
    err_white = geterrz_white(z, xx, cnt, dt, tottime); // approximation of err by the white-noise approximation
    err_res = initstd*initstd*pow(tottime, -2*z); // residual error
    // the real error should be err + err_res
    printf("%g\t%g\t%g\t%g\n", z, err, err_white, err_res);
  }
  return 0;
}
