/* adative velocity rescaling on a Lennard-Jones fluid */
#include "lj.h"
#include "ljeos.h"
#include <time.h>


int n = 108;
int nsteps = 500000;
double rho = 0.7;
double tp = 1.5;
double rcdef = 2.5;
double dt = 0.002;
const char *fnpos = "lj.pos";

int adaptive = 1; /* adaptive velocity scaling */
double fixene = -255.7;
double sfactor = 1.0; /* scaling factor */



/* print help message and die */
static void help(void)
{
  fprintf(stderr, "Adaptive velocity rescaling on a Lennard-Jones fluid\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -n:        set the number of particles, %d\n", n);
  fprintf(stderr, "  -E:        set the initial total energy, %g\n", fixene);
  fprintf(stderr, "  -S:        set the scaling factor, %g\n", sfactor);
  fprintf(stderr, "  -r:        set the density, %g\n", rho);
  fprintf(stderr, "  -T:        set the temperature, %g\n", tp);
  fprintf(stderr, "  -t:        set the number of steps, %d\n", nsteps);
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
      if ( strchr("nESrTtc", ch) != NULL ) {
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
          fixene = atof(q);
          adaptive = 0;
        } else if ( ch == 'S' ) {
          sfactor = atof(q);
        } else if ( ch == 'r' ) {
          rho = atof(q);
        } else if ( ch == 'T' ) {
          tp = atof(q);
        } else if ( ch == 't' ) {
          nsteps = atoi(q);
        } else if ( ch == 'c' ) {
          rcdef = atof(q);
        }
        break; /* skip the remaining characters */
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
void matchetot(lj_t *lj, int cycles, int stepspercycle)
{
  double ekin0 = 0, etot0, etot, etarget, de, s;
  int i, c, t;

  if ( adaptive ) {
    etarget = ljeos3d_get(rho, tp, NULL, NULL, NULL) * n
            + 0.5 * tp * lj->dof;
    /* add a random component */
    etarget += tp * sqrt(lj->dof) * randgaus();
  } else {
    etarget = fixene;
  }
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
  double bsm;
  double b2sm;
  double dbsm;
} betacm_t;


void accumbeta(lj_t *lj, betacm_t *acm)
{
  double bet = (lj->dof - 2) / (2 * lj->ekin);
  acm->cnt += 1;
  acm->bsm += bet;
  acm->b2sm += bet * bet;
  acm->dbsm += -bet / lj->ekin;
}


/* adaptively scale velocity */
void avscale(lj_t *lj, const betacm_t *acm,
    double *pdbde, double *pdbdk)
{
  double bet, dbet, dbdk, dbde, betav, betvar, de, s;
  int i;
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
  de = dbet / dbde;
  *pdbde = dbde;
  *pdbdk = dbdk;

  /* compute the scaling factor */
  s = sfactor * (de / lj->ekin) / acm->cnt;
  if ( s > 0.5 ) s = 0.5;
  else if ( s < -0.5 ) s = -0.5;
  s = sqrt(1 + s);
  /* scale the velocities */
  for ( i = 0; i < n; i++ )
    vsmul(lj->v[i], s);
}


int main(int argc, char **argv)
{
  int t;
  lj_t *lj;
  double cnt = 0, etot = 0, dbde = 0, dbdk = 0;
  double epsm = 0, etsm = 0, et2sm = 0;
  betacm_t acm[1] = {{0, 0, 0, 0}};

  doargs(argc, argv);

  mtscramble(time(NULL));
  lj = lj_open(n, rho, rcdef);
  lj->dof = n * 3 - 3;
  matchetot(lj, 3, 100);

  for ( t = 1; t <= nsteps; t++ ) {
    lj_vv(lj, dt);
    lj->ekin = lj_ekin(lj->v, n);
    if ( t % 1000 == 0 )
      printf("%8d\t%8.3f\t%8.3f\t%8.3f\n", t, lj->epot, lj->ekin, lj->epot + lj->ekin);
    if ( t > nsteps / 2 ) {
      cnt += 1;
      etot = lj->epot + lj->ekin;
      epsm += lj->epot;
      etsm += etot;
      et2sm += etot * etot;
    }
    accumbeta(lj, acm);
    /* adaptive velocity scaling */
    if ( adaptive ) avscale(lj, acm, &dbde, &dbdk);
  }
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);

  etot = lj->epot + lj->ekin;
  etsm /= cnt;
  et2sm = et2sm / cnt - etsm * etsm;
  fprintf(stderr, "rho %g, tp %g(%g), ep %g (ref. %g), "
      "etot %g, ave %g, var %g, dbde %g, dbdk %g, ratio %g\n",
      rho, tp, acm->cnt/acm->bsm, epsm/cnt/n,
      ljeos3d_get(rho, tp, NULL, NULL, NULL),
      etot, etsm, et2sm, dbde, dbdk, dbde/dbdk);
  return 0;
}

