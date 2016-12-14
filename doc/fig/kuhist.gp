#!/usr/bin/env gnuplot

# For the raw data, `make uhist` under data/wb

set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 7, 3.5 font "Times, 32"
set output "kuhist.eps"
set multiplot

rt = 0.543
ht = 0.88

set size rt, ht
set origin 0, 0
set rmargin 0

set xtics 100 offset 0, 0.3
set mxtics 10
#set xrange [-7720:-7400]
#set xlabel "Potential energy, {/Times-Italic U} (kcal/mol)" offset 0, 0
set xrange [1320:1540]
set xlabel "Kinetic energy, {/Times-Italic K} (kcal/mol)" offset 0, 0.5

set ytics 0.005
set mytics 5
set yrange [:0.018]
#set ylabel "Total energy, {/Times-Italic E}" offset 0.5, 0

#set key right Left reverse spacing 1.1 width -7
unset key

plot [:][:] \
    "../../data/wb/fix/Kfix.his"         u ($1):($2) w l  lt 1  lw 12.0 lc rgb "#c0c0c0" t "Microcanonical, {/Times-Italic NVE}", \
    "../../data/wb/adp/Kadp.his"         u ($1):($2) w l  lt 2  lw  3.0 t "Adaptive velocity scaling", \
    "../../data/wb/reg100/Kreg.his"      u ($1):($2) every 2 w p  pt 6  ps 1.5 lw  2.0 t "Regular velocity scaling", \
    "../../data/wb/can_vrs/Kcan_vrs.his" u ($1):($2) w l  lt 5  lw  5.0 t "Canonical, {/Times-Italic NVT}", \



set size 1 - rt, ht
set origin rt, 0
set lmargin 0.5
set rmargin 2.5

set xrange [-7720:-7430]
set xlabel "Potential energy, {/Times-Italic K} (kcal/mol)" offset 0, 0.5

set format y ""

unset key

plot [:][:] \
    "../../data/wb/fix/Ufix.his"         u ($1):($2) w l  lt 1  lw 12.0 lc rgb "#c0c0c0" t "Microcanonical, {/Times-Italic NVE}", \
    "../../data/wb/adp/Uadp.his"         u ($1):($2) w l  lt 2  lw  3.0 t "Adaptive velocity scaling", \
    "../../data/wb/reg100/Ureg.his"      u ($1):($2) every 2 w p  pt 6  ps 1.5 lw  2.0 t "Regular velocity scaling", \
    "../../data/wb/can_vrs/Ucan_vrs.his" u ($1):($2) w l  lt 5  lw  5.0 t "Canonical, {/Times-Italic NVT}", \


# dummy plots for the key
reset
set size 0.5, 1
set origin 0.06, ht - 0.05

set key below horizontal Left reverse maxrows 1 spacing 1.3
plot [] \
    x lt 1 lw 12.0 lc rgb "#c0c0c0" t "Microcanonical, {/Times-Italic NVE}", \
    x lt 2 lw  3.0 t "Adaptive velocity scaling"

set origin 0.51, ht - 0.05
plot [] \
    x w p pt 6 ps 2.0 lw  2.0 t "Regular velocity scaling", \
    x lt 5 lw  5.0 t "Canonical, {/Times-Italic NVT}"


unset multiplot
unset output
set terminal pop
reset
