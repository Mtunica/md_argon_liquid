set term png
reset

set linetype 1 lc rgb "black" lw 0.5 pt 0
set linetype 2 lc rgb "red" lw 0.5 pt 0
set linetype 3 lc rgb "blue" lw 0.5 pt 0
set linetype 4 lc rgb "green" lw 0.5 pt 0
set linetype 5 lc rgb "red" lw 2 pt 0
#-------------------------------------------
set xlabel("Time in reduced units")
set ylabel("Energy per atom in reduced units")

set output "data_results/energy_conservation.png"
set key inside top right

plot "data_results/estimation_simulation_1.dat" i 0 u 1:4 t "time step= 10^{-5}" with lp lt 1 ps 0,\
 	 "data_results/estimation_simulation_1.dat" i 1 u 1:4 t "time step=  10^{-4}" with lp lt 2 ps 0,\
 	 "data_results/estimation_simulation_1.dat" i 2 u 1:4 t "time step=  10^{-3}" with lp lt 3 ps 0

#-------------------------------------------
set xlabel("Time in reduced units")
set ylabel("Energy per atom in reduced units")

set output "data_results/energy_conservation2.png"
set key inside top left
plot "data_results/estimation_simulation_2.dat" i 0 u 1:4 t "time step= 10^{-5}" with lp lt 1 ps 0,\
 	 "data_results/estimation_simulation_2.dat" i 1 u 1:4 t "time step= 10^{-4}" with lp lt 2 ps 0
#-------------------------------------------

set xlabel("Time in reduced units")
set ylabel("Momentum per atom in reduced units")


set output "data_results/momentum_conservation.png"
set key inside top left
plot "data_results/estimation_simulation_1.dat" i 0 u 1:5 t "time step= 10^{-5}" with lp lt 1 ps 0,\
 	 "data_results/estimation_simulation_1.dat" i 1 u 1:5 t "time step= 10^{-4}" with lp lt 2 ps 0,\
 	 "data_results/estimation_simulation_1.dat" i 2 u 1:5 t "time step= 10^{-3}" with lp lt 3 ps 0
 
set output "data_results/momentum_conservation2.png"
set key inside top left
plot "data_results/estimation_simulation_2.dat" i 0 u 1:5 t "time step=10^{-6}" with lp lt 1 ps 0,\
 	 "data_results/estimation_simulation_2.dat" i 1 u 1:5 t "time step=10^{-5}" with lp lt 2 ps 0
 	 
#-------------------------------------------
set output "data_results/histogram0.png"
set style fill solid 1.0
set xlabel("Velocity in reduced units")
set ylabel ("Frequency over 1")
plot "data_results/histogram.dat" i 0 u 1:2 with boxes notitle

set output "data_results/histogram1.png"
set key inside top right
T=86.3

n= 966799239.0*sqrt(pi/2.0)*5*T*T*T/(39079286.0)
#n=1.0
f(x)=1.0/n*4.0*pi*(1.0/2.0*pi*T)**(3.0/2.0)*x*x*exp(-(x*x)/(2.0*T))


set style data histograms

plot "data_results/histogram.dat" i 1 u 1:2 with boxes t "Experimental data", f(x) lt 5 t "M-B distribution"

#-------------------------------------------

set xlabel("Time in reduced units")
set ylabel("Temperature in reduced units")

set output "data_results/temperature.png"

plot "data_results/temperature.dat" i 0 u 1:2 t "time step= 10^{-5}" with lp lt 1 ps 0,\
 	 "data_results/temperature.dat" i 1 u 1:2 t "time step=  10^{-4}" with lp lt 2 ps 0,\
 	 "data_results/temperature.dat" i 2 u 1:2 t "time step=  10^{-3}" with lp lt 3 ps 0


