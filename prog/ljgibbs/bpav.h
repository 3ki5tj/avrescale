typedef struct {
  double cnt, sbp, sdbp, sbvir, sbvir2;
  double sw, sw1, swbp;
  double bp, dbp, bvir, w, de, dvir, bp1;
  int n;
  double vol;
} bpav_t;


__inline void bpav_clear(bpav_t *bpav)
{
  bpav->cnt = 0;
  bpav->sbp = 0;
  bpav->sdbp = 0;
  bpav->sbvir = 0;
  bpav->sbvir2 = 0;
  bpav->sw = 0;
  bpav->swbp = 0;
}

/* get beta * pressure and the derivative d(beta*p)/d(1/vol) */
__inline static double lj_getbp(lj_t *lj, double beta,
    double *bvir, double *dbp)
{
  int n = lj->n;
  double bp, vol = lj->vol, rho, y;
  double irc, irc3, irc6, virtail, ytail;

  irc = 1/lj->rc;
  irc3 = irc * irc * irc;
  irc6 = irc3 * irc3;
  rho = n / vol;
  /* the tail correction for 3D */
  virtail = 32 * M_PI * rho * n / 3 * (irc6 - 1.5) * irc3;
  ytail = 32 * M_PI * rho * n * (4*irc6 - 3) * irc3;
  *bvir = beta * (lj->vir + virtail) / D;
  bp = (n + *bvir) / vol;
  y = 144 * lj->ep12 - 36 * lj->ep6 + ytail;
  *dbp = bp*vol + beta * y / (D*D);
  return bp;
}

/* energy and virial increase of adding a new particle */
__inline static double lj_denp(lj_t *lj, double *dvir)
{
  double xi[D], dx[D], l = lj->l, invl = 1/l, vol = lj->vol;
  double dr2, ir2, ir6, de, ep, ep6 = 0, ep12 = 0, eptail = 0;
  double irc, irc3, irc6, dvirtail, detail;
  int d, j, n = lj->n;

  /* randomly choose a position */
  for ( d = 0; d < D; d++ ) xi[d] = rand01() * lj->l;

  /* compute the energy change */
  for ( j = 0; j < n; j++ ) {
    dr2 = lj_pbcdist2(dx, xi, lj->x[j], l, invl);
    if ( dr2 >= lj->rc2 ) continue;
    ir2 = 1 / dr2;
    ir6 = ir2 * ir2 * ir2;
    ep12 += ir6 * ir6;
    ep6 += ir6;
  }
  ep6 *= 4;
  ep12 *= 4;
  ep = ep12 - ep6;
  *dvir = 12 * ep12 - 6 * ep6;

  /* compute the change of tail correction */
  eptail = lj_gettail(lj->rc, (n + 1)/lj->vol, n, NULL);
  de = ep + eptail - lj->epot_tail;

  /* add the tail correction */
  irc = 1 / lj->rc;
  irc3 = irc * irc * irc;
  irc6 = irc3 * irc3;
  dvirtail = (2 * n + 1) / vol * 32 * M_PI / 3 * (irc6 - 1.5) * irc3;
  detail = (2 * n + 1) / vol * 8 * M_PI / 9 * (irc6 - 3) * irc3;
  *dvir += dvirtail;
  de += detail;
  return de;
}

__inline void bpav_add(bpav_t *bpav, lj_t *lj, double beta)
{
  bpav->bp = lj_getbp(lj, beta, &bpav->bvir, &bpav->dbp);
  bpav->de = lj_denp(lj, &bpav->dvir);

  bpav->w = exp(-beta * bpav->de) * lj->vol / (lj->n + 1);
  bpav->bp1 = bpav->bp + bpav->dvir * beta / (D * lj->vol);

  if ( fabs(bpav->dbp) > 1e20 ) {
    fprintf(stderr, "dbp %g\n", bpav->dbp);
    exit(1);
  }

  bpav->cnt += 1;
  bpav->sbp += bpav->bp;
  bpav->sdbp += bpav->dbp;
  bpav->sbvir += bpav->bvir;
  bpav->sbvir2 += bpav->bvir * bpav->bvir;
  bpav->sw += bpav->w;
  bpav->swbp += bpav->w * bpav->bp1;

  /* save the number of particle and volume */
  bpav->n = lj->n;
  bpav->vol = lj->vol;
}

/* return the beta * pressure
 *   *dbp: d(beta*p)/d(rho)
 *   *w: < exp(-beta DU) >
 *   *dlnw: d(lnw)/d(lnV)
 * */
__inline double bpav_get(bpav_t *bpav, double *dbp, double *w, double *dlnw)
{
  int n = bpav->n;
  double cnt = bpav->cnt, bp, bvir, var, vol = bpav->vol;
  if ( cnt <= 0 ) {
    *dbp = 0;
    *w = 0;
    *dlnw = 0;
    return 0;
  }
  bp = bpav->sbp / cnt;
  bvir = bpav->sbvir / cnt;
  var = bpav->sbvir2 / cnt - bvir * bvir;
  *dbp = (bpav->sdbp / cnt - var) / n;
  *w = bpav->sw / cnt;
  *dlnw = 1 + (bpav->swbp / bpav->sw - bp) * vol;
  return bp;
}



