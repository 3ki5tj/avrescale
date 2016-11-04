#!/usr/bin/env python



'''
Compute the error versus the zooming factor
'''



import sys, os, getopt, shutil, re, math, glob
import zcom



fnout = None
cmdopt = ""
fnconf = ""
nsteps = None 
ntrials = 100
nproc = 1
zoom = 1.0
zrange = None
dKdE = 0
Etot = None
Edev = None
trace = False
thermostat = None
tNHCPeriod = 100.0
langRescaleDt = 20.0
verbose = 0



def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS]""" % sys.argv[0]

  print """
  Compute the error versus the zooming factor

  OPTIONS:

    -C, --conf=         set the template of the configuration file
    -t                  set the number of steps
    -M                  set the number of trials
    -p, --np=           set the number of processors
    -G, --dKdE=         set the explicit value of dKdE
    -Z, --zoom=         set the zooming factor, z
    -S, --scan=         set the range of the scanning z, format xmin:dx:xmax
    -E, --Etot=         set the target initial total energy
    -D, --Edev=         set the standard deviation of the initial total energy
    --trace             trace the time series of total energy
    --th=               set the thermostat method, Langevin, NH (Nose-Hoover), VRS (velocity rescaling)
    --tNHCPeriod        set the period of Nose-Hoover thermostat
    --langRescaleDt     set the time constant of the Langevin velocity rescaling thermostat
    --opt=              set options to be passed to the command line
    -o, --output=       set the output file
    -v                  be verbose
    --verbose=          set verbosity
    -h, --help          help
  """
  exit(1)



def doargs():
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "C:t:M:p:G:Z:S:E:D:hvo:",
        [
          "np=",
          "dKdE=", "zoom=", "scan=", "nsteps=",
          "cfg=", "conf=", "Etot=", "Edev=",
          "trace", "th=", "thermostat=",
          "tNHCPeriod=", "langRescaleDt=",
          "output=", "opt=",
          "help", "verbose=",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global fnconf, nsteps, ntrials, nproc, dKdE, zoom, zrange, Etot, Edev
  global trace, thermostat, tNHCPeriod, langRescaleDt
  global fnout, cmdopt, verbose

  for o, a in opts:
    if o in ("-C", "--cfg", "--conf",):
      fnconf = a
    elif o in ("-t", "--nsteps"):
      nsteps = int(a)
    elif o in ("-M",):
      ntrials = int(a)
    elif o in ("-p", "--np",):
      nproc = int(a)
    elif o in ("-G", "--dKdE",):
      dKdE = float(a)
    elif o in ("-E", "--Etot",):
      Etot = float(a)
    elif o in ("-D", "--Edev",):
      Edev = float(a)
    elif o in ("-Z", "--zoom",):
      zoom = float(a)
    elif o in ("-S", "--scan",):
      zrange = [float(x) for x in a.strip().split(":")]
    elif o in ("--trace",):
      trace = True
    elif o in ("--th", "--thermostat",):
      thermostat = a.upper()
    elif o in ("--tNHCPeriod",):
      tNHCPeriod = float(a)
    elif o in ("--langRescaleDt",):
      langRescaleDt = float(a)
    elif o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("--opt",):
      cmdopt = a
    elif o in ("-o", "--output"):
      fnout = a
    elif o in ("-h", "--help"):
      usage()


class Ave:
  ''' to compute average and variance '''
  def __init__(self):
    self.n = 0
    self.xsum = 0
    self.x2sum = 0

  def add(self, x):
    self.n += 1
    self.xsum += x
    self.x2sum += x * x

  def getave(self):
    return self.xsum / self.n

  def getvar(self):
    ave = self.xsum / self.n
    return self.x2sum / self.n - ave * ave



def getprogdir(build = True):
  ''' find the directory of NAMD '''
  progdir = "../../NAMD_mods/NAMD_2.11_thstat/Linux-x86_64-g++"
  i = 0
  while not os.path.isdir(progdir):
    progdir = "../" + progdir
    i += 1
    if i > 3: break
  if build: # build the program
    zcom.runcmd("make -C %s" % progdir)
  return progdir



def getfiles(initdir):
  ''' get files on the system '''
  global thermostat
  if not os.path.exists(initdir):
    initdir = "../" + initdir
  fnpsf = os.path.abspath( glob.glob(initdir + "/*.psf")[0] )
  fnpdb = fnpsf[:-3] + "pdb"
  fnprm = os.path.abspath( glob.glob(initdir + "/par_*")[0] )
  if fnconf:
    fncfg = fnconf
  else:
    # get 300Krand.conf or 300K.conf, not 300Kcan.conf or 300Kfix.conf
    ls = glob.glob(initdir + "/*Krand.conf")
    if len(ls):
      fncfg = ls[0]
    else:
      fncfg = glob.glob(initdir + "/*K.conf")[0]
  scfg = open(fncfg).readlines()
  for i in range(len(scfg)):
    ln = scfg[i].strip()
    if ln == "":
      scfg[i] = ""
    elif ( ln.startswith("rescaleAdaptive") or
           ln.startswith("langevin") or
           ln.startswith("energyLog") or
           ln.startswith("reinitvels") or
           ln.startswith("run") ):
      scfg[i] = ""
    elif ln.startswith("rescaleInitTotal") and Etot != None:
      scfg[i] = ""
    elif ln.startswith("rescaleInitDev") and Edev != None:
      scfg[i] = ""
    if thermostat:
      if ln.startswith("rescaleFreq") or ln.startswith("rescaleTemp"):
        scfg[i] = ""
  scfg = ''.join(scfg)
  return fnpsf, fnpdb, fnprm, scfg



def geterror(result):
  ''' get the initial and final errors from the output '''

  # construct the pattern
  numpat = "([\.0-9e+-]+)"
  pat = "etot %s, ave %s, var %s" % (numpat, numpat, numpat)
  m = re.search(pat, result)
  if not m:
    print "cannot find error information"
    exit(1)
  #print m.group(1), m.group(2), m.group(3)
  return float( m.group(1) )



def dosimul(zoom, build = True, fnlog = None):
  global fncfg, nsteps, cmdopt, Etot, Edev

  progdir = getprogdir(build)

  # copy files from the source directory
  fnpsf, fnpdb, fnprm, scfg = getfiles("normal")

  # create a temporary directory
  rundir = "tmprun"
  if not os.path.exists(rundir):
    os.mkdir(rundir)
  # move to the running directory
  os.chdir(rundir)

  # program directory
  prog = "../" + progdir + "/namd2 +p%s" % nproc

  # copy files
  os.system("cp %s ." % fnpsf)
  os.system("cp %s ." % fnpdb)
  os.system("cp %s ." % fnprm)

  if nsteps == None: nsteps = 10

  # command line
  fncfg = "run.conf"
  cmd = "%s %s %s" % (prog, cmdopt, fncfg)

  # output file
  if not fnlog: fnlog = "ez%s.log" % zoom
  fnlog = "../" + fnlog

  ln = "# zoom %s, nsteps %s, dKdE %s, Etot %s, Edev %s\n" % (
      zoom, nsteps, dKdE, Etot, Edev)
  open(fnlog, "a").write(ln)
  
  print "CMD: %s; LOG %s" % (cmd, fnlog)

  ef = Ave()
  ei = Ave()
  # multiple trials
  for i in range(ntrials):
    # write the configuration file
    fnene = "ene0.log"
    strcfg = scfg + '''
rescaleAdaptive           on
rescaleAdaptiveZoom       %s
rescaleAdaptiveFileFreq   1000
energyLogFile             %s
energyLogFreq             10
energyLogTotal            on
''' % (zoom, fnene)
    if dKdE > 0:
      strcfg += "rescaleAdaptiveDKdE       %s\n" % dKdE
    if Etot != None:
      strcfg += "rescaleInitTotal          %s\n" % Etot
    if Edev != None:
      strcfg += "rescaleInitDev            %s\n" % Edev

    strcfg += "run %s\n" % nsteps
    open(fncfg, "w").write(strcfg)

    # clear the directory
    os.system("rm -rf *.dat *.BAK *.old ene*.log *.vel *.xsc *.xst *.coor")

    ret, out, err = zcom.runcmd(cmd, capture = True, verbose = 0)
    # extract the energy
    ss = open(fnene).readlines()
    eitot = float( ss[0].split()[2] )
    j = -1
    while ss[j].strip() == "":
      j -= 1
    etot = float( ss[j].split()[2] )

    # update accumulators and print results
    ei.add(eitot)
    ef.add(etot)
    print "count %s, total energy %s, ave %s, var %s, init %s, iave %s, ivar %s" % (
        ef.n, etot, ef.getave(), ef.getvar(), eitot, ei.getave(), ei.getvar())

    ln = "%d %s %s\n" % (ef.n, etot, eitot)
    open(fnlog, "a").write(ln)

  # go back to the parent directory
  os.chdir("..")

  return ef.n, ef.getave(), ef.getvar(), ei.getave(), ei.getvar()



def dozscan():
  ''' scan the values of zoom factor '''
  global zrange
  zmin, zdel, zmax = zrange[0], zrange[1], zrange[2]
  sout = ""
  fnscan = "zscan.dat"
  zoom = zmin
  while zoom < zmax + zdel * 0.5:
    cnt, eave, evar, eiave, eivar = dosimul(zoom, zoom == zmin)
    ln = "%s %s %s %s %s %s\n" % (
        zoom, eave, evar, eiave, eivar, cnt)
    sout += ln
    print ln,
    open(fnscan, "w").write(sout)
    zoom += zdel
  print "saved results to", fnscan



def dotrace(zoom, build = True, fntrace = None):
  global fncfg, nsteps, cmdopt, Etot, Edev
  global thermostat, tNHCPeriod, langRescaleDt

  progdir = getprogdir(build)

  if Edev == None: Edev = 100
  if Etot == None: Etot = -6141.16
  Etot += Edev
  Edev = 0

  # copy files from the source directory
  fnpsf, fnpdb, fnprm, scfg = getfiles("normal")

  # create a temporary directory
  rundir = "tmprun"
  if not os.path.exists(rundir):
    os.mkdir(rundir)
  # move to the running directory
  os.chdir(rundir)

  # program directory
  prog = "../" + progdir + "/namd2 +p%s" % nproc

  # copy files
  os.system("cp %s ." % fnpsf)
  os.system("cp %s ." % fnpdb)
  os.system("cp %s ." % fnprm)

  if nsteps == None: nsteps = 1

  # command line
  fncfg = "run.conf"
  cmd = "%s %s %s" % (prog, cmdopt, fncfg)

  # output file
  if not fntrace: fntrace = "etot.tr"
  fntrace = "../" + fntrace

  print "CMD: %s; TRACE %s" % (cmd, fntrace)

  eparr = []
  ekarr = []
  etarr = []
  # multiple trials
  for i in range(ntrials):
    # write the configuration file
    fnene = "ene0.log"

    strcfg = scfg + '''
energyLogFile             %s
energyLogFreq             1
energyLogTotal            on
''' % (fnene)

    if not thermostat:
      strcfg += '''
rescaleAdaptive           on
rescaleAdaptiveZoom       %s
rescaleAdaptiveFileFreq   1000
''' % (zoom)
      if dKdE > 0:
        strcfg += "rescaleAdaptiveDKdE       %s\n" % dKdE

    elif thermostat == "NH":
      # Nose-Hoover thermostat
      strcfg += '''
tNHC              on
tNHCTemp          $temperature
tNHCLen           5
tNHCPeriod        %g
''' % (tNHCPeriod)

    elif thermostat.startswith("VR"):
      # Langevin-style velocity rescaling
      strcfg += '''
langRescale       on
langRescaleTemp   $temperature
langRescaleDt     %g
''' % (langRescaleDt)

    else: # if thermostat.startswith("L"):
      # Langevin thermostat
      strcfg += '''
langevin          on
langevinDamping   1
langevinTemp      $temperature
'''

    if Etot != None:
      strcfg += "rescaleInitTotal          %s\n" % Etot
    if Edev != None:
      strcfg += "rescaleInitDev            %s\n" % Edev

    strcfg += "run %s\n" % nsteps
    open(fncfg, "w").write(strcfg)

    # clear the directory
    os.system("rm -rf *.dat *.BAK *.old ene*.log *.vel *.xsc *.xst *.coor")

    ret, out, err = zcom.runcmd(cmd, capture = True, verbose = 0)
    # extract the energy
    ss = open(fnene).read().strip().split('\n')
    ssplit = [ln.strip().split() for ln in ss]
    tm = len(ss)

    # update accumulators and print results
    if i == 0:
      etarr = [ Ave() for k in range(tm) ]
      eparr = [ Ave() for k in range(tm) ]
      ekarr = [ Ave() for k in range(tm) ]
    for k in range(tm):
      epot = float( ssplit[k][1] )
      eparr[k].add(epot)
      etot = float( ssplit[k][2] )
      etarr[k].add(etot)
      ekin = etot - epot
      ekarr[k].add(ekin)

    print "count %s, total energy, final: ave %s, fvar %s, init: ave %s, var %s" % (
        i + 1, etarr[-1].getave(), etarr[-1].getvar(),
        etarr[0].getave(), etarr[0].getvar())

    strace = "# count %s, zoom %s, nsteps %s, dKdE %s, Etot %s, Edev %s\n" % (
      etarr[0].n, zoom, nsteps, dKdE, Etot, Edev)
    for k in range(tm):
      strace += "%s\t%g\t%g\t%g\t%g\t%g\t%g\n" % ( ssplit[k][0],
          etarr[k].getave(), etarr[k].getvar(),
          eparr[k].getave(), eparr[k].getvar(),
          ekarr[k].getave(), ekarr[k].getvar() )
    open(fntrace, "w").write(strace)

  # go back to the parent directory
  os.chdir("..")



if __name__ == "__main__":
  doargs()
  if trace:
    dotrace(zoom, fntrace = fnout)
  elif zrange:
    dozscan()
  else:
    dosimul(zoom, fnlog = fnout)

