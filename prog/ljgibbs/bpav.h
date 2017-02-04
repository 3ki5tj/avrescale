typedef struct {
  double cnt, sbpV, sdbp, sbvir, sbvir2;
  double sw, swbpV;
  double bp, dbp, bvir, w, de, dvir, bpV1;
} bpav_t;


__inline void bpav_clear(bpav_t *bpav)
{
  bpav->cnt = 0;
  bpav->sbpV = 0;
  bpav->sdbp = 0;
  bpav->sbvir = 0;
  bpav->sbvir2 = 0;
  bpav->sw = 0;
  bpav->swbpV = 0;
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

__inline static void bpav_add(bpav_t *bpav, lj_t *lj, double beta)
{
  double vol = lj->vol;

  bpav->bp = lj_getbp(lj, beta, &bpav->bvir, &bpav->dbp);
  bpav->de = lj_denp(lj, &bpav->dvir);

  bpav->w = exp(-beta * bpav->de) * vol / (lj->n + 1);
  bpav->bpV1 = bpav->bp * vol + beta * bpav->dvir / D;

  if ( fabs(bpav->dbp) > 1e20 ) {
    fprintf(stderr, "dbp %g\n", bpav->dbp);
    exit(1);
  }

  bpav->cnt += 1;
  bpav->sbpV += bpav->bp * vol;
  bpav->sdbp += bpav->dbp;
  bpav->sbvir += bpav->bvir;
  bpav->sbvir2 += bpav->bvir * bpav->bvir;
  bpav->sw += bpav->w;
  bpav->swbpV += bpav->w * bpav->bpV1;
}

/* return the beta * pressure
 *   *dbp: d(beta*p)/d(rho)
 *   *w: < exp(-beta DU) >
 *   *dlnw: d(lnw)/d(lnV)
 * */
__inline double bpav_get(bpav_t *bpav, lj_t *lj, double *dbp, double *w, double *dlnw)
{
  double cnt = bpav->cnt, bpV, bvir, var, vol = lj->vol;

  if ( cnt <= 0 ) {
    *dbp = 0;
    *w = 0;
    *dlnw = 0;
    return 0;
  }
  bpV = bpav->sbpV / cnt;
  bvir = bpav->sbvir / cnt;
  var = bpav->sbvir2 / cnt - bvir * bvir;
  *dbp = (bpav->sdbp / cnt - var) / lj->n;
  *w = bpav->sw / cnt;
  *dlnw = 1 + bpav->swbpV / bpav->sw - bpV;
  return bpV / vol;
}



