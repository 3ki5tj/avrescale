# T: number of steps
# a: initial error (as the variance of the energy)
# g: Gamma, correlation integral
f(z,T,a,g)=a/T**(2*z)+g*z*z/(2*z-1)*(T**(2*z-1)-1)/T**(2*z)

# to be read from the first line of the output of program md
# the var value, or specified by the option -Ddev, aval = dev*dev
aval = 722.25

# To be obtained from constant a magnitude run of the program md
# Example 1:
#   ./md -A 0.001 -G0.72 -t 100000000 -E-255.7 -D1e-100 > mag0.001.log
#   ../mkhist.py --col=4 mag0.001.log --t0=1000000 --dx=0.0001
# if the var is 3.2, then Gamma = 2*3.2/0.001 = 6400
# Example 2:
#   ./md -A 0.0005 -G0.72 -t 200000000 -E-255.7 -D1e-100 > mag0.0005.log
#   ../mkhist.py --col=4 mag0.001.log --t0=1000000 --dx=0.0001
# if the var is 1.6, then Gamma = 2*1.6/0.0005 = 6400
gval = 6400

#set logscale y

plot [0.4:] [3:30]\
  f(x,1000,aval,gval), "lj1/zscan.dat" u 1:3, \
  f(x,2000,aval,gval), "lj3/zscan.dat" u 1:3

# To show the relative accuracy
#plot [] \
#  "lj1/zscan.dat" u 1:($3/f($1,1000,aval,gval)) w lp, \
#  "lj3/zscan.dat" u 1:($3/f($1,2000,aval,gval)) w lp, 1
