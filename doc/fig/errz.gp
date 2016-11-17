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
# with rescaleAdaptiveMag 0.0001
# namd 300Kmag.conf
# ../mkhist.py ene0mag.log --col=3 --dx=0.002
# If the value of 'var' is 18
# Then Gamma = 2*18/0.0001 = 360000
gval = 273000

set logscale x
#set xtics 1 offset 0, 0.3
#set mxtics 2
#set xrange [0.05:10]
set xlabel "{/Times-Italic z}" offset 0, 0.5

set logscale y
set format y "10^{/*0.8 %T}"
set yrange [0.1:400]
#set ytics 5
#set mytics 5
set ylabel "Error of the total energy" offset 1, 0

# `width` to reduce the text length
set key right Left reverse spacing 1.2 width -3

plot [0.25:][:] \
    "../../data/wb/wb_t100000.dat"  u ($1):2  w p  pt 7  ps 2.0 t "Simulation, {/Times-Italic T} = 10^{/*0.8 5}", \
    f(x,100000,aval,gval)                     lt 1 lw 3         t "Prediction, {/Times-Italic T} = 10^{/*0.8 5}", \
    "../../data/wb/wb_t1M.dat"      u ($1):2  w p  pt 5  ps 1.8 t "Simulation, {/Times-Italic T} = 10^{/*0.8 6}", \
    f(x,1000000,aval,gval)                    lt 3 lw 3         t "Prediction, {/Times-Italic T} = 10^{/*0.8 6}"


unset multiplot
unset output
set terminal pop
reset
