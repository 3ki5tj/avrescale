#!/usr/bin/env gnuplot

set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 5, 3.5 font "Times, 32"
set output "errz.eps"
set multiplot



# error vs the proportionality constant z
# of the inverse-time schedule of the updating magnitude
#   alpha(t) = z / t

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
# If the value of var is 136.5
# Then Gamma = 2*136.5/0.001 = 273000
gval = 273000

set logscale x
#set xtics 1 offset 0, 0.3
#set mxtics 2
#set xrange [0.05:10]
set xlabel "{/Times-Italic z}" offset 0, 0.5

#set logscale y
#set format y "10^{/*0.8 %T}"
set ytics 5
set mytics 5
#set yrange [5e-7:0.1]
set ylabel "Error" offset 1, 0

# `width` to reduce the text length
set key right Left reverse spacing 1.2

plot [0.3:][:20] \
    "../../data/wb_t100000.dat" u ($1):2         w p  pt 7  ps 2.0 t "Simulation", \
    f(x,100000,aval,gval) lt 1 t "Prediction"


unset multiplot
unset output
set terminal pop
reset
