#!/usr/bin/env gnuplot

set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 5, 3.5 font "Times, 32"
set output "uhist.eps"
set multiplot

set xtics 100 offset 0, 0.3
set mxtics 5
set xrange [-7720:-7400]
set xlabel "Potential energy, {/Times-Italic U} (kcal/mol)" offset 0, 0

set ytics 0.01
set mytics 5
set yrange [:0.026]
#set ylabel "Total energy, {/Times-Italic E}" offset 0.5, 0

set key right Left reverse spacing 1.1 width -7

plot [:][:] \
    "../../data/waterbox_can/ene0can.his" u ($1):($2) w l  lt 1  lw 2.0 t "Canonical", \
    "../../data/waterbox_fix/ene0fix.his" u ($1):($2) w l  lt 2  lw 2.0 t "Microcanonical", \
    "../../data/waterbox_adp/ene0adp.his" u ($1):($2) w l  lt 4  lw 2.0 t "Adaptive velocity scaling", \

# "../../data/waterbox_reg/ene0reg.his" u ($1):($2) w l  lt 5  lw 2.0 t "Regular velocity scaling"

unset multiplot
unset output
set terminal pop
reset
