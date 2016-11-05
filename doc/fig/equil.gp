#!/usr/bin/env gnuplot

set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 5, 3.7 font "Times, 32"
set output "equil.eps"
set multiplot

set xtics 100 offset 0, 0.3
set mxtics 2
set xlabel "Number of MD steps, {/Times-Italic t}" offset 0, 0.3

set ytics 20
set mytics 2
#set yrange [:-5980]
set ylabel "Total energy, {/Times-Italic E}" offset 1.5, 0

set key right Left reverse spacing 1.1 width -8

plot [:200][:] \
    "../../data/wb/trace/set1/etot.tr"          u ($1):($6) w lp  lt 1 pt 7  ps 1.0 t "Adaptive velocity scaling", \
    "../../data/wb/trace/set1/etot_Lang.tr"     u ($1):($6) w lp  lt 2 pt 4  ps 1.0 t "Langevin dynamics", \
    "../../data/wb/trace/set1/etot_NHP100.tr"   u ($1):($6) w lp  lt 5 pt 8  ps 1.2 t "Nos\351-Hoover chain", \
    "../../data/wb/trace/set1/etot_VRSDt30.tr"  u ($1):($6) w lp  lt 6 pt 10 ps 1.2 t "Velocity rescaling", \


unset multiplot
unset output
set terminal pop
reset
