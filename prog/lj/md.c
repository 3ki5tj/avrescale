/* adative velocity rescaling on a Lennard-Jones fluid */
#include "lj.h"
#include "ljeos.h"
#include <time.h>


int n = 108;
long nsteps = 500000;

double rho = 0.7;
double tp = 1.5;
double rcdef = 2.5;
double dt = 0.002;
const char *fnpos = "lj.pos";
int dopr = 0;

int adaptive = 1; /* adaptive velocity scaling */
int surf = 0; /* use the surface temperature */
double initene = 0; /* -255.7; */
double initdev = 0;
double zfactor = 1.0; /* zooming factor */
double mag = 0.0;
double dKdE = 0.0;

int nstlog = 1000; /* logging frequency */


/* print help message and die */
static void help(void)
{
  fprintf(stderr, "Adaptive velocity rescaling on a Lennard-Jones fluid\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -n:        set the number of particles, %d\n", n);
  fprintf(stderr, "  -F:        turn off adaptive averaging, %d\n", !adaptive);
  fprintf(stderr, "  -S:        use the surface temperature, %d\n", surf);
  fprintf(stderr, "  -E:        set the initial total energy, %g\n", initene);
  fprintf(stderr, "  -D:        set the standard deviation of the initial total energy, %g\n", initdev);
  fprintf(stderr, "  -Z:        set the zooming factor for the scaling magnitude, %g\n", zfactor);
  fprintf(stderr, "  -G:        set the explicit value of dKdE, %g\n", dKdE);
  fprintf(stderr, "  -A:        set the fixed scaling magnitude, %g\n", mag);
  fprintf(stderr, "  -g:        set the logging interval, %d\n", nstlog);
  fprintf(stderr, "\n === [ Lennard-Jones Parameters ] ===\n\n");
  fprintf(stderr, "  -r:        set the density, %g\n", rho);
  fprintf(stderr, "  -T:        set the temperature, %g\n", tp);
  fprintf(stderr, "  -t:        set the number of steps, %ld\n", nsteps);
  fprintf(stderr, "  -c:        set the cutoff distance, %g\n", rcdef);
  fprintf(stderr, "  -h:        display this message\n");
  exit(1);
}

/* parse command-line arguments */
static int doargs(int argc, char **argv)
{
  int i, j, ch;
  char *p, *q;

  for ( i = 1; i < argc; i++ ) {
    /* test if the token is an argument or an option */
    if ( argv[i][0] != '-' ) continue;

    /* test long argument, like --help */
    if ( argv[i][1] == '-' ) {
      p = argv[i] + 2;
      /* let q point to the argument of the option */
      if ( (q = strchr(p, '=')) != NULL ) {
        *q++ = '\0';
      } else {
        q = NULL;
      }
      /* do matching */
      if ( strcmp(p, "help") == 0 ) {
        help();
      } else {
        fprintf(stderr, "Unknown long option, key [%s], val [%s]\n",
            p, (q != NULL ? q : "NULL") );
      }
      continue;
    }

    /* this is a short option */
    for ( j = 1; (ch = argv[i][j]) != '\0'; j++ ) {
      if ( strchr("nEDZAGrTtcg", ch) != NULL ) {
        q = p = argv[i] + j + 1;
        if ( *p != '\0' ) {
          /* the argument follows the option immediately
           * e.g., -n100 */
          q = p;
        } else if ( ++i < argc ) {
          /* the argument is the next token of command line */
          q = argv[i];
        } else {
          fprintf(stderr, "Option -%c requires an argument!\n", ch);
          help();
        }

        if ( ch == 'n' ) {
          n = atoi(q);
        } else if ( ch == 'E' ) {
          initene = atof(q);
        } else if ( ch == 'D' ) {
          initdev = atof(q);
        } else if ( ch == 'Z' ) {
          zfactor = atof(q);
        } else if ( ch == 'G' ) {
          dKdE = atof(q);
        } else if ( ch == 'A' ) {
          mag = atof(q);
        } else if ( ch == 'r' ) {
          rho = atof(q);
        } else if ( ch == 'T' ) {
          tp = atof(q);
        } else if ( ch == 't' ) {
          nsteps = atol(q);
        } else if ( ch == 'c' ) {
          rcdef = atof(q);
        } else if ( ch == 'g' ) {
          nstlog = atoi(q);
        }
        break; /* skip the remaining characters */
      } else if ( ch == 'F' ) {
        adaptive = 0;
      } else if ( ch == 'S' ) {
        surf = 1;
      } else if ( ch == 'h' ) {
        help();
      } else {
        fprintf(stderr, "Unknown short option %s, j %d, ch %c\n",
            argv[i], j, ch);
        help();
      }
    }
  }
  return 0;
}


/* scale the velocities to reach the desired total energy */
static void matchetot(lj_t *lj, int cycles, int stepspercycle)
{
  double ekin0 = 0, etot0, etot, etarget, de, s, std = 0;
  int i, c, t;

  if ( initene != 0 ) {
    etarget = initene;
  } else {
    etarget = ljeos3d_get(rho, tp, NULL, NULL, NULL) * n
            + 0.5 * tp * lj->dof;
  }
  if ( initdev != 0 ) {
    std = initdev;
  } else {
    std = tp * sqrt(lj->dof);
  }
  /* add a random component */
  etarget += std * randgaus();

  fprintf(stderr, "matching %s energy %g, std %g, var %g; ",
      (adaptive ? "random" : "fixed"), etarget, std, std * std);
  fprintf(stderr, "n %d, dof %d, rho %g, T %g\n", n, lj->dof, rho, tp);
  for ( c = 1; c <= cycles; c++ ) {
    /* compute the average total energy over a few steps */
    etot0 = 0;
    for ( t = 0; t < stepspercycle; t++ ) {
      lj_vv(lj, dt);
      ekin0 = lj_ekin(lj->v, n);
      etot0 += lj->epot + ekin0;
    }
    etot0 /= stepspercycle;
    de = etarget - etot0;
    s = sqrt(1 + de / ekin0);
    /* scale the velocities */
    for ( i = 0; i < n; i++ )
      vsmul(lj->v[i], s);
    lj->ekin = lj_ekin(lj->v, n);
    etot = lj->epot + lj->ekin;
    fprintf(stderr, "cycle %d, epot %8.3f, ekin %8.3f -> %8.3f, etot %8.3f\n",
        c, lj->epot, ekin0, lj->ekin, etot);
  }
}



/* accumulator of the temperature data */
typedef struct {
  double cnt;
  double tpsm;
  double tp2sm;
  double bsm;
  double b2sm;
  double dbsm;
} betacm_t;


static void accumbeta(lj_t *lj, betacm_t *acm)
{
  double bet = (lj->dof - 2) / (2 * lj->ekin);
  double tp1 = 2 * lj->ekin / lj->dof;
  acm->cnt += 1;
  acm->tpsm += tp1;
  acm->tp2sm += tp1 * tp1;
  acm->bsm += bet;
  acm->b2sm += bet * bet;
  acm->dbsm += -bet / lj->ekin;
}

static double getkin(lj_t *lj, betacm_t *acm, double *var)
{
  double tpave, tpvar, c;
  tpave = acm->tpsm / acm->cnt;
  tpvar = acm->tp2sm / acm->cnt - tpave * tpave;
  c = lj->dof / 2;
  *var = tpvar * c * c;
  return tpave * c;
}

/* adaptively scale velocity */
static void avscale(lj_t *lj, const betacm_t *acm, int surface,
    double *pdbde, double *pdbdk, double *pdkde)
{
  double kin, dkin, tpav, betav, dkde, dkde0;
  double bet, dbet, dbdk, dbde, betvar;
  double de, s;
  int i;

  if ( surface ) { /* surface temperature */
    bet = (lj->dof - 2) / (2 * lj->ekin);
    dbet = 1/tp - bet;
    betav = acm->bsm / acm->cnt;
    betvar = acm->b2sm / acm->cnt - betav * betav;
    dbdk = acm->dbsm / acm->cnt;
    if ( acm->cnt < 10 ) {
      dbde = dbdk;
    } else {
      dbde = dbdk + betvar;
      if ( dbde > 0.1 * dbdk ) dbde = 0.1 * dbdk;
    }
    if ( dKdE > 0 ) {
      dbde = -dKdE/(0.5*lj->dof + dKdE - 2)/(tp*tp);
    }
    de = dbet / dbde;
    *pdbde = dbde;
    *pdbdk = dbdk;
    *pdkde = dbde / dbdk;

  } else { /* bulk temperature */

    kin = 0.5 * lj->dof * tp;
    dkin = kin - lj->ekin;
    tpav = acm->tpsm / acm->cnt;
    betav = acm->bsm / acm->cnt;
    dkde0 = (1 - betav * tpav) * 0.5 * lj->dof;
    if ( dKdE > 0 ) {
      dkde = dKdE;
    } else if ( acm->cnt < 10 ) {
      dkde = 1;
    } else {
      dkde = dkde0;
      if ( dkde < 0.1 ) dkde = 0.1;
    }
    de = dkin / dkde;
    *pdbde = -0.5 * lj->dof / kin / kin * dkde0;
    *pdkde = dkde0;
    *pdbdk = *pdbde / *pdkde;
  }

  /* compute the scaling factor */
  if ( mag > 0 ) {
    s = mag * (de / lj->ekin);
  } else {
    s = zfactor * (de / lj->ekin) / acm->cnt;
  }
  if ( s > 0.5 ) s = 0.5;
  else if ( s < -0.5 ) s = -0.5;
  s = sqrt(1 + s);
  /* scale the velocities */
  for ( i = 0; i < n; i++ )
    vsmul(lj->v[i], s);
}


int main(int argc, char **argv)
{
  long t;
  lj_t *lj;
  double cnt = 0, etot = 0, dbde = 0, dbdk = 0, dkde = 0;
  double epsm = 0, etsm = 0, et2sm = 0, eitot;
  double epref, kinave, kinvar;
  betacm_t acm[1] = {{0}};

  doargs(argc, argv);

  mtscramble(time(NULL) + 314159 * clock());
  lj = lj_open(n, rho, rcdef, dopr);
  lj->dof = n * 3 - 3;
  matchetot(lj, 3, 100);
  eitot = lj->epot + lj->ekin;
  printf("%8d\t%8.3f\t%8.3f\t%8.3f\n", 0, lj->epot, lj->ekin, lj->epot + lj->ekin);

  for ( t = 1; t <= nsteps; t++ ) {
    lj_vv(lj, dt);
    lj->ekin = lj_ekin(lj->v, n);
    if ( t % nstlog == 0 )
      printf("%8ld\t%8.3f\t%8.3f\t%8.3f\n", t, lj->epot, lj->ekin, lj->epot + lj->ekin);
    if ( t > nsteps / 2 ) {
      cnt += 1;
      etot = lj->epot + lj->ekin;
      epsm += lj->epot;
      etsm += etot;
      et2sm += etot * etot;
    }
    accumbeta(lj, acm);
    /* adaptive velocity scaling */
    if ( adaptive ) avscale(lj, acm, surf, &dbde, &dbdk, &dkde);
  }

  etot = lj->epot + lj->ekin;
  etsm /= cnt;
  et2sm = et2sm / cnt - etsm * etsm;
  epref = ljeos3d_get(rho, tp, NULL, NULL, NULL);
  kinave = getkin(lj, acm, &kinvar);
  fprintf(stderr, "rho %g, tp %g(%g,%g), ep %g (ref. %g), "
      "etot %g, ave %g, var %g, einit %g, kin %g, var %g",
      rho, tp, acm->tpsm/acm->cnt, acm->cnt/acm->bsm, epsm/cnt/n, epref,
      etot, etsm, et2sm, eitot, kinave, kinvar);
  fprintf(stderr, ", db/dE %g, db/dK %g, dK/dE %g", dbde, dbdk, dkde);
  if ( mag > 0 ) fprintf(stderr, ", Gamma %g", 2*et2sm/mag);
  fprintf(stderr, "\n");

  lj_close(lj);
  return 0;
}

