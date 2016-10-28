#!/usr/bin/env python



'''
'''



import sys, os, getopt, shutil, re, math, glob
import zcom



fnout = None
cmdopt = ""
fnconf = ""
nsteps = 10000
ntrials = 100
nproc = 4
zoom = 1.0
zrange = None
dKdE = 0
Etot = None
Edev = None
verbose = 0



def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS]""" % sys.argv[0]

  print """
  Compute the error over a range of zooming factor

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
    --opt=              set options to be passed to the command line
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
          "output=", "opt=",
          "help", "verbose=",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global fnconf, nsteps, ntrials, nproc, dKdE, zoom, zrange, Etot, Edev
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



def getfiles(initdir):
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



def dosimul(zoom, build = True):
  global fncfg, fnout, cmdopt, Etot, Edev

  # build the program
  progdir = "../../NAMD_mods/NAMD_2.11_thstat/Linux-x86_64-g++"
  if not os.path.isdir(progdir):
    progdir = "../" + progdir
  if build:
    zcom.runcmd("make -C %s" % progdir)

  # copy files from the source directory
  fnpsf, fnpdb, fnprm, scfg = getfiles("waterbox")

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

  # command line
  fncfg = "run.conf"
  cmd = "%s %s %s" % (prog, cmdopt, fncfg)

  fnlog = "../ez%s.log" % zoom
  ln = "# zoom %s, nsteps %s, dKdE %s, Etot %s, Edev %s\n" % (
      zoom, nsteps, dKdE, Etot, Edev)
  open(fnlog, "a").write(ln)
  
  print "CMD: %s; LOG %s" % (cmd, fnlog)

  cnt = 0
  esum = 0
  e2sum = 0
  eisum = 0
  ei2sum = 0
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

    # print results
    cnt += 1
    eisum += eitot
    ei2sum += eitot * eitot
    eiave = eisum / cnt
    eivar = ei2sum / cnt - eiave * eiave
    esum += etot
    e2sum += etot * etot
    eave = esum / cnt
    evar = e2sum / cnt - eave * eave
    print "count %s, total energy %s, ave %s, var %s, init %s, iave %s, ivar %s" % (
        cnt, etot, eave, evar, eitot, eiave, eivar)

    ln = "%d %s %s\n" % (cnt, etot, eitot)
    open(fnlog, "a").write(ln)

  # go back to the parent directory
  os.chdir("..")

  return cnt, eave, evar, eiave, eivar



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



if __name__ == "__main__":
  doargs()
  if zrange:
    dozscan()
  else:
    dosimul(zoom)

