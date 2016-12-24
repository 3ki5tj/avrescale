typedef struct {
  double cnt, sumx, sumx2;
} av_t;

static void av_add(av_t *av, double x)
{
  av->cnt += 1;
  av->sumx += x;
  av->sumx2 += x * x;
}

static double av_getave(av_t *av)
{
  return av->cnt > 0 ? av->sumx / av->cnt : 0;
}

static void av_clear(av_t *av)
{
  av->cnt = 0;
  av->sumx = 0;
  av->sumx2 = 0;
}

#define av_init(av) av_clear(av)
