f(z,T,a,g)=a/T**(2*z)+g*z*z/(2*z-1)*(T**(2*z-1)-1)/T**(2*z)

# to be read from the first line of the output of program md
# the var value
aval = 722.25

# number of steps
Tval = 2000

# estimated value
gval = 6000

plot f(x,Tval,aval,gval), "lj3/zscan2.dat" u 1:3
