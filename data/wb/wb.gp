# T: number of steps
# a: initial error (as the variance of the energy)
# g: Gamma, correlation integral
f(z,T,a,g)=a/T**(2*z)+g*z*z/(2*z-1)*(T**(2*z-1)-1)/T**(2*z)

# To be specified by the option in the configuration file
# rescaleInitDev 100
aval = 100000

# To be obtained from constant a magnitude simulation, 300Kmag.conf
# with rescaleAdaptiveMag 0.0001
# namd 300Kmag.conf
# ../mkhist.py ene0mag.log --col=3 --dx=0.002
# If the value of var is 18
# Then Gamma = 2*18/0.0001 = 360000
gval = 360000

set logscale

plot [0.3:][:60] \
  f(x,1000000,aval,gval) t "Prediction, T = 10^6", \
  "wberrz_t1M.dat" u 1:2 w l t "Prediction, T = 10^6", \
  "wb_t1M.dat" u 1:($2) t "Simulation, T = 10^6", \
  f(x,100000,aval,gval) t "Prediction, T = 10^5", \
  "wberrz_t100000.dat" u 1:2 w l t "Prediction, T = 10^5", \
  "wb_t100000.dat" u 1:($2) t "Simulation, T = 10^5"



# wb_t100000 are created by
#../runnamd.py -p1 -t100000 --dKdE=0.33 -Z0.3 -M10000
#../runnamd.py -p1 -t100000 --dKdE=0.33 -Z0.5 -M10000
#../runnamd.py -p1 -t100000 --dKdE=0.33 -Z1 -M10000
#../runnamd.py -p1 -t100000 --dKdE=0.33 -Z3 -M10000
#../runnamd.py -p1 -t100000 --dKdE=0.33 -Z10 -M10000
# After run
#   mv run*/*.log wb_t100000
# For the final error
# ../mkhist.py --col=2 --dx=0.001
# For the initial error
# ../mkhist.py --col=3 --dx=0.001
