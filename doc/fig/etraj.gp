#!/usr/bin/env gnuplot

set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced size 5, 3.7 font "Times, 32"
set output "etraj.eps"
set multiplot

#set xtics 1 offset 0, 0.3
#set mxtics 2
set xlabel "Number of MD steps, {/Times-Italic t} ({/Symbol \264} 10^{/*0.7 5 })" offset 0, 0.3

set ytics 10
set mytics 10
set yrange [:-6100]
set ylabel "Total energy, {/Times-Italic E}" offset 0.3, 0

set key right Left reverse spacing 1.2 width -7

plot [:30][:] \
    "../../data/wb/fix/ene0fix.log" u ($1/1e5):($3) every 5000 w p  lt 1 pt 6  ps 1.0 t "Without velocity scaling", \
    "../../data/wb/reg/ene0reg.log" u ($1/1e5):($3) every 5000 w p  lt 1 pt 1  ps 1.0 t "Regular velocity scaling", \
    "../../data/wb/adp/ene0adp.log" u ($1/1e5):($3) every 5000 w p  lt 2 pt 9  ps 1.0 t "Adaptive velocity scaling", \

#plot [:5][:] \
#    "../../data/waterbox_fix/set1/ene0fix.log" u ($1/1e5):($3) every 500  w p  pt 1  ps 1.0 t "Regular MD", \
#    "../../data/waterbox_adp/set1/ene0adp.log" u ($1/1e5):($3) every 500  w p  pt 7  ps 1.0 t "Adaptive velocity scaling", \

unset multiplot
unset output
set terminal pop
reset
