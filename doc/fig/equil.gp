#!/usr/bin/env gnuplot

set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 5, 3.7 font "Times, 32"
set output "equil.eps"
set multiplot

set xtics 200 offset 0, 0.3
set mxtics 2
set xlabel "Number of MD steps, {/Times-Italic t}" offset 0, 0.5

set ytics 50 offset 0.3, 0
set mytics 5
#set yrange [:-5980]
set ylabel "Total energy, {/Times-Italic E}" offset 1.2, 0

set key right Left reverse spacing 1.1 width -8

plot [:1000][:-6000] \
    "../../data/wb/trace/set2/etot_AVS.tr"      u ($1):($2) every 20 w  p  pt 6  ps 1.2 notitle, \
    "../../data/wb/trace/set2/etot_Lang.tr"     u ($1):($2) every 20 w  p  pt 4  ps 1.0 notitle, \
    "../../data/wb/trace/set2/etot_NHCP100.tr"  u ($1):($2) every 20 w  p  pt 8  ps 1.2 notitle, \
    "../../data/wb/trace/set2/etot_VRSDt100.tr" u ($1):($2) every 20 w  p  pt 10 ps 1.2 notitle, \
    "../../data/wb/trace/set2/etot_AVS.tr"      u ($1):($2) w l   lt 1              notitle, \
    "../../data/wb/trace/set2/etot_Lang.tr"     u ($1):($2) w l   lt 2              notitle, \
    "../../data/wb/trace/set2/etot_NHCP100.tr"  u ($1):($2) w l   lt 5              notitle, \
    "../../data/wb/trace/set2/etot_VRSDt100.tr" u ($1):($2) w l   lt 6              notitle, \
    0                                                       w lp  lt 1 pt 7  ps 1.0 t "Adaptive velocity scaling", \
    0                                                       w lp  lt 2 pt 4  ps 1.0 t "Langevin dynamics", \
    0                                                       w lp  lt 5 pt 8  ps 1.2 t "Nos\351-Hoover chain", \
    0                                                       w lp  lt 6 pt 10 ps 1.2 t "Velocity rescaling"


unset multiplot
unset output
set terminal pop
reset
