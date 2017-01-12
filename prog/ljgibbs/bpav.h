typedef struct {
  double cnt, sbp, sdbp, sbvir, sbvir2;
} bpav_t;


void bpav_clear(bpav_t *bpav)
{
  bpav->cnt = 0;
  bpav->sbp = 0;
  bpav->sdbp = 0;
  bpav->sbvir = 0;
  bpav->sbvir2 = 0;
}

void bpav_add(bpav_t *bpav, lj_t *lj, double beta)
{
  int n = lj->n;
  double bp, bvir, dbp, vol = lj->vol, rho, y;
  double irc, irc3, irc6, virtail, ytail;

  irc = 1/lj->rc;
  irc3 = irc * irc * irc;
  irc6 = irc3 * irc3;
  rho = n / vol;
  /* the tail correction for 3D */
  virtail = 32 * M_PI * rho * n/3 * (irc6 - 1.5) * irc3;
  ytail = 32 * M_PI * rho * n * (4*irc6 - 3) * irc3;
  bvir = beta * (lj->vir + virtail) / (D * vol);
  bp = n / vol + bvir;
  y = 144 * lj->ep12 - 36 * lj->ep6 + ytail;
  dbp = -bp/vol - beta * y / (D*D*vol*vol);

  bpav->cnt += 1;
  bpav->sbp += bp;
  bpav->sdbp += dbp;
  bpav->sbvir += bvir;
  bpav->sbvir2 += bvir * bvir;
}

double bpav_get(bpav_t *bpav, double *dbp)
{
  double cnt = bpav->cnt, bp, bvir, var;
  if ( cnt <= 0 ) {
    *dbp = 0;
    return 0;
  }
  bp = bpav->sbp / cnt;
  bvir = bpav->sbvir / cnt;
  var = bpav->sbvir2 / cnt - bvir * bvir;
  *dbp = bpav->sdbp / cnt + var;
  return bp;
}



