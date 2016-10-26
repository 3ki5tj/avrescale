#!/usr/bin/env python



'''
'''



import sys, os, getopt, shutil, re, math
import zcom



fnout = None
cmdopt = ""
nsteps = 10000
ntrials = 100
zoom = 1.0
zrange = None
dKdE = 0
verbose = 0



def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS]""" % sys.argv[0]

  print """
  Compute the error over a range of the zooming factor, z

  OPTIONS:

    -t                  set the number of steps
    -M                  set the number of trials
    -G, --dKdE=         set the explicit value of dKdE
    -Z, --zoom=         set the zooming factor, z
    -S, --scan=         set the range of scanning z, format xmin:dx:xmax
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
        "t:M:G:Z:S:hvo:",
        [
          "dKdE=", "zoom=", "scan=",
          "output=", "opt=",
          "help", "verbose=",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global nsteps, ntrials, dKdE, zoom, zrange
  global fnout, cmdopt, verbose

  for o, a in opts:
    if o in ("-t",):
      nsteps = int(a)
    elif o in ("-M",):
      ntrials = int(a)
    elif o in ("-G", "--dKdE",):
      dKdE = float(a)
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
  return float( m.group(1) ), float( m.group(2) )



def dosimul(zoom, build = True):
  global fncfg, fnout, cmdopt

  # build the program
  progdir = "../prog/lj"
  if not os.path.isdir(progdir):
    progdir = "../" + progdir
  if build:
    zcom.runcmd("make -C %s" % progdir)

  prog = "md"
  try:
    # make a copy of the program in case
    # it gets modified or deleted later
    shutil.copy("%s/%s" % (progdir, prog), "./%s" % prog)
  except:
    pass

  cmd = "./%s -t%s -Z%s %s" % (prog, nsteps, zoom, cmdopt)
  if dKdE > 0: cmd += "-G%s" % dKdE
  cmd = cmd.strip()

  fnlog = "ez%s.log" % zoom
  ln = "# zoom %s, nsteps %s, dKdE %s\n" % (
      zoom, nsteps, dKdE)
  open(fnlog, "w").write(ln) # empty the log
  
  print "CMD: %s; LOG %s" % (cmd, fnlog)

  cnt = 0
  esum = 0
  e2sum = 0
  for i in range(ntrials):
    ret, out, err = zcom.runcmd(cmd, capture = True, verbose = 0)
    # get the last line of the output
    result = err.strip().split('\n')[-1]
    etot, etav = geterror(result)

    # print results
    cnt += 1
    esum += etot
    e2sum += etot * etot
    eave = esum / cnt
    evar = e2sum / cnt - eave * eave
    print "count %s, total energy %s, ave %s, var %s" % (
        cnt, etot, eave, evar)

    ln = "%d %s %s\n" % (cnt, etot, etav)
    open(fnlog, "a").write(ln)

  return cnt, eave, evar


def dozscan():
  global zrange
  zmin, zdel, zmax = zrange[0], zrange[1], zrange[2]
  sout = "# nsteps %d dKdE %s %s\n" % (
      nsteps, dKdE, cmdopt)
  fnscan = "zscan.dat"
  zoom = zmin
  while zoom < zmax + zdel * 0.5:
    cnt, eave, evar = dosimul(zoom, zoom == zmin)
    ln = "%s %s %s %s\n" % (zoom, eave, evar, cnt)
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
