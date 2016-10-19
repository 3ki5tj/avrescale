#!/usr/bin/env python



'''
'''



import sys, os, getopt, shutil, re, math, glob
import zcom



fnout = None
cmdopt = ""
nsteps = 10000
ntrials = 100
nproc = 4
eta = 1.0
etarange = None
dKdE = 0
verbose = 0



def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS]""" % sys.argv[0]

  print """
  Compute the error over a range of eta 

  OPTIONS:

    -t                  set the number of steps
    -M                  set the number of trials
    -p, --np=           set the number of processors
    --dKdE=             set the explicit value of dKdE
    --eta=              set the multiplier
    -S, --scan=         set the range of the scanning eta, format xmin:dx:xmax
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
        "t:M:p:S:hvo:",
        [
          "np=",
          "dKdE=", "eta=", "scan=", "nsteps=",
          "output=", "opt=",
          "help", "verbose=",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global nsteps, ntrials, nproc, dKdE, eta, etarange
  global fnout, cmdopt, verbose

  for o, a in opts:
    if o in ("-t", "--nsteps="):
      nsteps = int(a)
    elif o in ("-M",):
      ntrials = int(a)
    elif o in ("-p", "--np",):
      nproc = int(a)
    elif o in ("--dKdE",):
      dKdE = float(a)
    elif o in ("--eta",):
      eta = float(a)
    elif o in ("-S", "--scan",):
      etarange = [float(x) for x in a.strip().split(":")]
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
  # get 300K.conf not 300Kcan.conf or 300Kfix.conf
  fncfg = glob.glob(initdir + "/*K.conf")[0]
  scfg = open(fncfg).readlines()
  for i in range(len(scfg)):
    if ( scfg[i].startswith("rescaleAdaptive") or
         scfg[i].startswith("langevin") or
         scfg[i].startswith("energyLog") or
         scfg[i].startswith("run") ):
      scfg[i] = "\n"
    if scfg[i].strip() == "":
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



def dosimul(eta, build = True):
  global fncfg, fnout, cmdopt 

  if build:
    # build the program
    progdir = "../../NAMD_mods/NAMD_2.11_thstat/Linux-x86_64-g++"
    if not os.path.isdir(progdir):
      progdir = "../" + progdir
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

  cnt = 0
  esum = 0
  e2sum = 0
  # multiple trials
  for i in range(ntrials):
    # clear the directory
    os.system("rm -rf *.dat *.BAK *.old *.log *.vel *.xsc *.xst *.coor")

    # write the configuration file
    fncfg = "run.conf"
    fnene = "ene0.log"
    strcfg = scfg + '''
rescaleAdaptive           on
rescaleAdaptiveEta        %s
rescaleAdaptiveFileFreq   1000
energyLogFile             %s
energyLogFreq             10
energyLogTotal            on
''' % (eta, fnene)
    if dKdE > 0:
      strcfg += "rescaleAdaptiveDKdE %s\n" % dKdE
    strcfg += "run %s\n" % nsteps
    open(fncfg, "w").write(strcfg)

    # command line
    cmd = "%s %s %s" % (prog, cmdopt, fncfg)
    if i == 0: print cmd

    ret, out, err = zcom.runcmd(cmd, capture = True, verbose = 0)
    # extract the energy
    arr = open(fnene).readlines()[-1].split()
    etot = float(arr[2])

    # print results
    cnt += 1
    esum += etot
    e2sum += etot * etot
    eave = esum / cnt
    evar = e2sum / cnt - eave * eave
    print "count %s, total energy %s, ave %s, var %s" % (
        cnt, etot, eave, evar)

  # go back to the parent directory
  os.chdir("..")

  return cnt, eave, evar



def doscan():
  ''' scan the values of eta '''
  global etarange
  etamin = etarange[0]
  etadel = etarange[1]
  etamax = etarange[2]
  eta = etamin
  sout = ""
  while eta < etamax + etadel * 0.5:
    cnt, eave, evar = dosimul(eta, eta == etamin)
    ln = "%s %s %s %s\n" % (eta, eave, evar, cnt)
    sout += ln
    print ln,
    eta += etadel
  open("scan.out", "w").write(sout)



if __name__ == "__main__":
  doargs()
  if etarange:
    doscan()
  else:
    dosimul()

