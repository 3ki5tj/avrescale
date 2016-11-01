#!/usr/bin/env gnuplot

set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 5, 3.5 font "Times, 32"
set output "etraj.eps"
set multiplot

#set xtics 1 offset 0, 0.3
#set mxtics 2
set xlabel "Number of MD steps, {/Times-Italic t} ({/Symbol \264} 10^{/*0.7 5 })" offset 0, 0

set ytics 5
set mytics 5
set yrange [:-6135]
set ylabel "Total energy, {/Times-Italic E}" offset 0.5, 0

set key right Left reverse spacing 1.2 width -5

plot [:10][:] \
    "../../data/waterbox_fix/ene0fix.log" u ($1/1e5):($3) every 500  w p  pt 1  ps 1.0 t "Regular MD", \
    "../../data/waterbox_adp/ene0adp.log" u ($1/1e5):($3) every 500  w p  pt 7  ps 1.0 t "Adaptive velocity scaling", \

#plot [:5][:] \
#    "../../data/waterbox_fix/set1/ene0fix.log" u ($1/1e5):($3) every 500  w p  pt 1  ps 1.0 t "Regular MD", \
#    "../../data/waterbox_adp/set1/ene0adp.log" u ($1/1e5):($3) every 500  w p  pt 7  ps 1.0 t "Adaptive velocity scaling", \

unset multiplot
unset output
set terminal pop
reset
