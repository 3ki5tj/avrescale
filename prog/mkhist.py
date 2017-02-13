#!/usr/bin/env python



import sys, os, glob, getopt
from math import *



fnins = ["ene.log",]
fnout = None
dx = 0.5
dE = 0.5
dT = 1
T = 300
col = 2
colT = -1
colE = -1
fnrst = ""
t0 = 0
i0 = 0
every = 1
tpinterp = 10
drop1 = False
docorr = False
fncorr = None
corrmax = 200
nohist = False


def showhelp():
  print "mkhist.py [Options] ene.log"
  print "Description:"
  print "  Make a histogram of an energy log file"
  print "  It can compute the autocorrelation function with --corr"
  print "Options:"
  print "  -T, --tp=:           set the temperature"
  print "  -d, --dx=:           set the energy bin size"
  print "  --dT=:               set the temperature tolerance"
  print "  -o, --output=:       set the output file"
  print "  -i, --input=:        set the input file"
  print "  -c, --col=:          set the column for the quantity, can also be --col=dif for the difference between column 3 and 2, or --col=invK for the inverse (for autocorrelations)"
  print "  --t0=:               set the first time step"
  print "  --i0=:               set the first frame"
  print "  -e, --every=:        set the stride of frames"
  print "  --colT=:             set the column for temperature"
  print "  --rst=:              set the restart file for WHAM (weighted histogram analysis method)"
  print "  --tm=:               set the number of subdivisions for each temperature bin (for WHAM)"
  print "  --colE=:             set the column for energy (for reweighting)"
  print "  --dE=:               set the bin size for the energy grid (for reweighting)"
  print "  -1, --drop1          drop the last frame"
  print "  --corr, --acf        compute the autocorrelation function"
  print "  --fncorr=, --fnacf=  set the output file name for the autocorrelation function"
  print "  --corrmax=           set the maximal number of steps for the autocorrelation function"
  print "  --nohist             skipping histogram"
  exit(1)



def doargs():
  global fnins, fnout, dx, dE, dT, T, col, colT, colE, fnrst
  global t0, i0, every, tpinterp, drop1
  global docorr, fncorr, corrmax, nohist

  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "hd:T:o:i:c:e:1",
        ["dx=", "dE=", "de=", "dT=", "dt=", "tp=", "input=", "output=",
         "col=", "colT=", "colE=", "rst=", "t0=", "tm=", "drop=",
         "corr", "acf", "fncorr=", "fnacf=", "corrmax=", "nohist"])
  except:
    print "Error parsing the command line"
    sys.exit(1)

  fnins = args
  for o, a in opts:
    if o == "-h":
      showhelp()
    elif o in ("-T", "--tp"):
      T = float(a)
    elif o in ("--dT", "--dt"):
      dT = float(a)
    elif o in ("-d", "--dx"):
      dx = float(a)
    elif o in ("--dE", "--de"):
      dE = float(a)
    elif o in ("-i", "--input"):
      fnins += [a,]
    elif o in ("-o", "--output"):
      fnout = a
    elif o in ("-c", "--col"):
      col = a.lower() # keep the string form, for "dif" or "invk"
      if col.isdigit():
        col = int(col)
    elif o in ("--colT",):
      colT = int(a)
    elif o in ("--colE",):
      colE = int(a)
    elif o in ("--rst",):
      fnrst = a
    elif o in ("--t0",):
      t0 = float(a)
    elif o in ("--i0",):
      i0 = int(a)
    elif o in ("-e", "--every",):
      every = int(a)
    elif o in ("--tm",):
      tm = int(a)
    elif o in ("-1", "--drop1",):
      drop1 = True
    elif o in ("--corr", "--acf",):
      docorr = True
    elif o in ("--fncorr", "--fnacf",):
      fncorr = a
      docorr = True
    elif o in ("--corrmax",):
      corrmax = int(a)
      docorr = True
    elif o in ("--nohist",):
      nohist = True

  if len(fnins) == 0:
    fnins = glob.glob("e*.log")
  #print args, opts, fnins, fnout, T, dT, dE
  if not fnins: showhelp()



class Hist:
  def __init__(self, xmin1, xmax1, dx):
    if xmin1 >= xmax1:
      xmax1 = xmin1 + 10 * dx
      xmin1 -= 10 * dx
    self.dx = dx
    self.n = int( (xmax1 - xmin1) / dx + 0.9999 )
    self.xmin = round(xmin1 / dx) * dx
    self.xmax = self.xmin + self.n * self.dx
    self.arr = [0] * self.n
    self.sum1 = 1e-300
    self.sumx = 0
    self.sumx2 = 0

  def add(self, newx, w = 1.0):
    if newx < self.xmin:
      # extend the histogram array to the left
      lftx = newx - self.dx * 10 # leave some margin
      deln = int( (self.xmin - lftx) / self.dx )
      delx = deln * self.dx
      self.xmin -= delx
      self.n += deln
      self.arr = [0] * deln + self.arr

    if newx >= self.xmax - self.dx:
      # extend the histogram array to the right
      rtx = newx + self.dx * 10
      deln = int( (rtx - self.xmax) / self.dx )
      delx = deln * self.dx
      self.xmax += delx
      self.n += deln
      self.arr += [0] * deln

    i = int( (newx - self.xmin) / self.dx + 0.5 )
    if i < 0 or i >= self.n:
      print i, self.n, newx, self.xmin, self.xmax
      raw_input()
    self.arr[i] += w

    self.sum1 += w
    self.sumx += w * newx
    self.sumx2 += w * newx * newx

  def save(self, fn):
    fac = 1.0 / (self.dx * self.sum1)
    # determine the boundaries
    i0 = 0
    while i0 < self.n:
      if self.arr[i0] > 0: break
      i0 += 1
    i1 = self.n - 1
    while i1 >= i0:
      if self.arr[i1] > 0: break
      i1 -= 1

    i = i0
    s = ""
    sumx = 0
    sumx2 = 0
    while i <= i1:
      dist = self.arr[i] * fac
      x = self.xmin + (i + 0.5) * self.dx
      sumx += x * self.arr[i]
      sumx2 += x * x * self.arr[i]
      s += "%s\t%s\t%s\n" % (x, dist, self.arr[i])
      i += 1
    avx = self.sumx / self.sum1
    varx = self.sumx2 / self.sum1 - avx * avx
    print "saving histogram file %s, total %s, ave %s, var %s, std %s" % (
        fn, self.sum1, avx, varx, sqrt(varx) )
    if fn.lower() not in ("null", "none"):
      open(fn, "w").write(s)



# Boltzmann constant in kcal/mol in NAMD
kB = 0.001987191

def lnadd(x, y):
  ''' return ln(exp(x) + exp(y)) = x + ln(1 + exp(y-x)) '''
  if x < y: # swap
    x, y = y, x
  return x + log(1 + exp(y -x))


class WHAM:
  def __init__(self, Tref, fn):
    ''' construct the WHAM object from a NAMD restart file '''
    self.bref = 1/(kB * Tref)
    try:
      s = open(fn).readlines()
      self.is_open = True
    except:
      if fn: print "cannot open %s, will use explicit filter with dT = %s" % (fn, dT)
      self.is_open = False
      return
    info = s[0].split()
    self.bmin = 1.0 / (kB * float(info[2]))
    self.bmax = 1.0 / (kB * float(info[3]))
    ntp0 = int(info[4])
    dbeta0 = (self.bmax - self.bmin) / ntp0
    self.eav0 = [0] * ntp0
    self.lnw0 = [0] * ntp0
    xp = 0
    for i in range(ntp0):
      ln = s[i + 1].split()
      self.eav0[i] = float(ln[1])
      invw = float(ln[7])
      xp += invw / (self.bmin + dbeta0 *(i+.5))
    # assuming the ensemble weight is beta^(-xp)
    xp = round(xp / ntp0, 3)
    # integrate to get lnz
    self.lnz0 = [0] * (ntp0 + 1)
    for i in range(ntp0):
      self.lnz0[i + 1] = self.lnz0[i] - self.eav0[i] * dbeta0
    self.interp(ntp0, tpinterp, xp)
    self.lnwarr = None

  def interp(self, ntp0, m, xp):
    ''' compute the partition function and weights of the finer grids '''
    self.ntp = ntp0 * m
    self.dbeta = (self.bmax - self.bmin) / self.ntp
    self.bmid = [0] * self.ntp
    self.lnwb = [0] * self.ntp
    self.lnzmid = [0] * self.ntp
    self.lnz = [0] * (self.ntp + 1)
    for i in range(self.ntp):
      self.bmid[i] = self.bmin + self.dbeta * (i + 0.5)
      self.lnwb[i] = -xp * log(self.bmid[i])
      j = i % m
      k = i / m
      x = (j + 0.5) / m
      self.lnzmid[i] = self.lnz0[k] * (1 - x) + self.lnz0[k + 1] * x
      x = (j + 1.0) / m
      self.lnz[i + 1] = self.lnz0[k] * (1 - x) + self.lnz0[k + 1] * x
      #print i, self.bmid[i], self.lnwb[i], self.lnz[i], self.lnzmid[i]
    #raw_input()

  def getlnw(self, bet, ene):
    ''' directly compute the weight for energy at temperature bet '''
    # compute the numerator
    x = (bet - self.bmin) / self.dbeta
    ib = int( x )
    dx = x - ib
    lnzb = self.lnz[ib]*(1 - dx) + self.lnz[ib + 1]*dx
    lnnum = -bet * ene - lnzb
    # compute the denominator
    lnden = -100000
    for i in range(self.ntp):
      xp = -self.bmid[i] * ene + self.lnwb[i] - self.lnzmid[i]
      lnden = lnadd(lnden, xp)
    lnden += log(self.dbeta)
    return lnnum - lnden

  def getlnwarr(self, emin, emax, de):
    ''' tabulate the logarithm of weight for the energy grid '''
    self.de = de
    imin = int(emin / de) - 1
    imax = int(emax / de) + 1
    self.en = imax - imin
    self.emin = imin * de
    self.emax = imax * de
    self.lnwarr = [0] * (self.en + 1)
    print "tabulating WHAM weight from %s to %s, %d bins ... " % (self.emin, self.emax, self.en),
    for i in range(self.en + 1):
      ene = self.emin + i * de
      self.lnwarr[i] = self.getlnw(self.bref, ene)
      #print ene, self.lnwarr[i]
    print "done!"
    #raw_input("%s %s %s" % (emin, emax, de))

  def getweight(self, ene):
    ''' get the ensemble weight for the given energy '''
    if not self.lnwarr: # initialize
      emin = ene - 100 * dE
      emax = ene + 100 * dE
      self.getlnwarr(emin, emax, dE)
    elif ene < self.emin: # extend to the left
      emin = ene - 100 * dE
      self.getlnwarr(emin, self.emax, dE)
    elif ene > self.emax: # extend to the right
      emax = ene + 100 * dE
      self.getlnwarr(self.emin, emax, dE)

    # linearly interpolate the weight
    x = (ene - self.emin) / self.de
    i = int( x )
    dx = x - i
    lnwt = self.lnwarr[i] * (1 - dx) + self.lnwarr[i + 1] * dx
    return exp(lnwt)



def getcol(col, arr):
  ''' interpret the column data '''
  if type(col) == str:
    if col.startswith("dif") or col in ("k", "kin", "kinetic"):
      return arr[2] - arr[1]
    elif col in ("ik", "invk"):
      return 1.0/(arr[2] - arr[1])
    else:
      print "unknown column type: '%s'" % col
      raise Exception
  else:
    return arr[col-1]



def mkhist_simple(s, fnout):
  global col, i0, every
  n = len(s)
  if drop1: n -= 1 # drop the last frame
  hist = None
  for i in range(n):
    ln = s[i].strip()
    # skip a comment line
    if ln == "" or ln.startswith("#"): continue
    i += 1
    if i <= i0 or i % every != 0: continue
    tok = [float(y) for y in ln.split()]
    try:
      tm = tok[0]
    except Exception:
      print "error in %s: %s" % (fn, ln)
      break
    if tm <= t0: continue
    x = getcol(col, tok)
    if not hist:
      hist = Hist(x, x, dx)
    hist.add(x)
  if hist:
    hist.save(fnout)



def mkhist_reweight(s, fnout):
  ''' temperature selected/reweighted histogram '''
  global col, colT, colE, fnrst, i0, every
  n = len(s) - 1 # drop the last frame
  hist = None
  wham = WHAM(T, fnrst) # make a WHAM object
  for i in range(n):
    ln = s[i].strip()
    # skip a comment line
    if ln == "" or ln.startswith("#"): continue
    i += 1
    if i <= i0 or i % every != 0: continue
    tok = [float(x) for x in ln.split()]
    try:
      tm = tok[0]
      tp = tok[colT - 1]
    except:
      break
    if tm <= t0: continue
    x = getcol(col, tok)
    ene = float(tok[colE - 1])
    if wham.is_open: # use WHAM
      w = wham.getweight(ene)
    else: # use temperature-bin filter
      if tp < T - dT or tp > T + dT: continue
      w = 1.0
    if not hist:
      hist = Hist(x, x, dx)
    hist.add(x, w)
  if hist:
    hist.save(fnout)



def mkhist(fnin):
  global fnout, colE, colT

  s = open(fnin).readlines()
  fno = fnout
  if not fno:
    fno = os.path.splitext(fnin)[0] + ".his"

  if col == 2:
    # automatically determine file type
    for ln in s:
      if ln.startswith("#"): continue
      arr = ln.split()
      num = len( arr )
      if num == 2:
        colT = 0
        colE = 2
      else:
        temp = float(arr[2])
        if num == 3 and (temp >= 100 and temp <= 600):
          # for a WHAM simulation the format should be 
          #  step energy temperature
          # so the third column should be between 100K to 600K
          colT = 3
          colE = 2
        elif num == 4: # step x temperature energy
          colT = 3
          colE = 4
      break

  if colT > 0:
    mkhist_reweight(s, fno)
  else:
    mkhist_simple(s, fno)



def getcorr(fnin):
  global col, corrmax, t0

  t = []
  x = []
  i = 0
  for ln0 in open(fnin).readlines():
    ln = ln0.strip()
    # skip a comment line
    if ln == "" or ln.startswith("#"): continue
    i += 1
    if i <= i0 or i % every != 0: continue
    tok = [float(y) for y in ln.strip().split()]
    tm = int(tok[0])
    if tm <= t0: continue
    xv = getcol(col, tok)
    x.append( xv )
    if len(t) < 2: t.append( tm )
  dt = t[1] - t[0]

  n = len(x)
  var = 0
  ave = sum(x)/n
  dx = [y - ave for y in x]
  xx = []
  imax = min(corrmax + 1, n - 1)
  for i in range(imax):
    # compute autocoorelation of i steps apart
    xxv = sum(dx[j]*dx[j+i] for j in range(n - i))/(n - i)
    xx.append(xxv)
  var = xx[0]

  # compute the autocorrelation integral
  ss = 0
  sout = ""
  good = True
  for i in range(len(xx) - 1):
    y1 = xx[i]
    y2 = xx[i+1]
    if y1 < 0: good = False
    if good: ss += 0.5*(y2+y1)
    # assuming exponential variation within two neighboring points
    #ss += (y2-y1)/log(y2/y1)
    sout += "%s %s\n" % (i*dt, xx[i]/xx[0])
  ss *= 2*dt
  # autocorrelation time
  tau = (ss / var - 1) / 2

  # add a tag line
  sout = "# %s %s %s %s %s %s %s\n" % (n, ave, var, dt, len(xx) - 1, ss, tau) + sout

  if type(col) == str:
    scol = col
  else:
    scol = "col" + str(col)

  if fncorr != None:
    fnout = fncorr
  else:
    fnout = os.path.splitext(fnin)[0] + "_" + scol + ".acf"

  open(fnout, "w").write(sout)
  print "saved correlation function to %s, n %s, ave %s, var %s, " % (
      fnout, n, ave, var),
  print "correlation integral %s, autocorrelation time %s" % (
      ss, tau)



if __name__ == "__main__":
  doargs()
  for fn in fnins:
    if not nohist: mkhist(fn)
    if docorr: getcorr(fn)
