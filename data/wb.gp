# T: number of steps
# a: initial error (as the variance of the energy)
# g: Gamma, correlation integral
f(z,T,a,g)=a/T**(2*z)+g*z*z/(2*z-1)*(T**(2*z-1)-1)/T**(2*z)

# To be specified by the option in the configuration file
# rescaleInitDev 100
aval = 100000

# To be obtained from constant a magnitude simulation, 300Kmag.conf
# with rescaleAdaptiveMag 0.001
# namd 300Kmag.conf
# ../mkhist.py ene0mag.log --col=3 --dx=0.002
# If the value of var is 137
# Then Gamma = 2*137/0.001 = 274000
gval = 274000

#set logscale y

plot [0.3:][:200] \
  f(x,5000,aval,gval), "wb_t5000.dat" u 1:($2)

# wb_t5000 are created by
#../runnamd.py -p2 -t5000 --dKdE=0.33 -Z1 -M10000
#../runnamd.py -p1 -t5000 --dKdE=0.33 -Z0.3 -M10000
#../runnamd.py -p2 -t5000 --dKdE=0.33 -Z3 -M10000
#../runnamd.py -p1 -t5000 --dKdE=0.33 -Z4 -M10000
#../runnamd.py -p2 -t5000 --dKdE=0.33 -Z0.5 -M10000
# After run
#   mv run*/*.log wb_t5000
# For the final error
# ../mkhist.py --col=2 --dx=0.001
# For the initial error
# ../mkhist.py --col=3 --dx=0.001
