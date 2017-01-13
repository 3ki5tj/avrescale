# load the equation of state
load "ljeosKN.gp"
#load "ljeosJZG.gp"

reset

tp = 1.15
beta = 1 / tp
pres = 0.061

Fr(rho, T) = Fex(rho, T) + T*(log(rho) + 1)
mu(rho, T) = muex(rho, T) + T + T*(log(rho) + 1)

set samples 500, 500
#plot [1.0:40] Fex(1/x, tp) - 0.1*x
#plot [1.6:20] P(1/x, tp), pres
#plot [1.6:20] beta * Fex(1/x, tp)
#plot [0.05:0.7] Fr(x, tp) + pres/x, P(x, tp) axis x1y2, pres axis x1y2

#Gibbs Ensemble
n = 200
vol = 700
Ftot(x, y) = n*x*Fr(n*x/(vol*y), tp) + n*(1-x)*Fr(n*(1-x)/(vol*(1-y)), tp)
Fcut = -n*1.4
Ftrunc(x, y) = (Ftot(x, y) > Fcut) ? Fcut : Ftot(x, y)
#plot [0.08:0.92] Ftrunc(0.5, x)

set xlabel "N"
set ylabel "V"
#set isosamples 20, 20
#splot [0.01:0.99][0.08:0.92] Ftrunc(x, y)
set pm3d map
set isosamples 100, 100
splot [0.0:n][0.08*vol:0.92*vol] Ftrunc(x/n, y/vol) w pm3d

#n = 100
#vol = 700
#F2(v1, v2) = Fr(n/v1, tp) + Fr(n/v2, tp)
#F2max = -2.5
#F2trunc(v1, v2) = F2(v1, v2) > F2max ? F2max : F2(v1, v2)
#set xlabel "V1"
#set ylabel "V2"
#splot [0.1:0.9][0.1:0.9] F2trunc(vol*x, vol*y)

