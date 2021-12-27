reset

set term png


set linetype 1 lc rgb "black" lw 0.5 pt 0
set linetype 2 lc rgb "red" lw 0.5 pt 0
set linetype 3 lc rgb "blue" lw 0.5 pt 0
set linetype 4 lc rgb "green" lw 0.5 pt 0

set style line 5 lc rgb 'black' pt 7 ps 1.5 lt 1 lw 2 
set style line 6 lc rgb 'red' pt 7 ps 1.5 lt 1 lw 2 
set style line 7 lc rgb 'blue' pt 7 ps 1.5 lt 1 lw 2 


set output "data_results/displacement.png"
set ylabel("MSD (Angstroms^2)")
set xlabel("Time (ps)")

set key inside top left

g(x) = 6.0*D*x + b


fit g(x) "data_results/mean_square_displacement.dat" u ($1/10):2 via D, b


plot "data_results/mean_square_displacement.dat" u ($1/10):2 notitle with lp lt 1 ps 2, g(x) t "MSD=6D*t"



set output "data_results/g.png"

set linetype 9 lc rgb "black" lw 0.5 pt 0.5

set ylabel("Radial function")
set xlabel("Radius (Angstroms)")

plot "data_results/histogram_b.dat" u 1:2 notitle with lp lt 9 ps 0.5


set output "data_results/energies.png"

set key inside bottom left

set xlabel("Density (g/cm^3)")
set ylabel("Energy per atom (Kj/mol)")

g(x)=A*x+B
h(x)=C*x+D
t(x)=E*x+F


fit g(x)  "data_results/enegies.dat" u 1:2 via A,B
fit h(x)  "data_results/enegies.dat" u 1:3 via C,D
fit t(x)  "data_results/enegies.dat" u 1:4 via E,F

plot "data_results/enegies.dat" u 1:2 t "Kinetic energy" with points lc rgb 'black' pt 7 ps 1, g(x) with lp lt 1 ps 0 notitle,\
	 "data_results/enegies.dat" u 1:3 t "Potential energy"  with points lc rgb 'red' pt 7 ps 1, h(x)  with lp lt 2 ps 0 notitle,\
	 "data_results/enegies.dat" u 1:4 t "Total energy"  with points lc rgb 'blue' pt 7 ps 1, t(x)  with lp lt 3 ps 0 notitle
	 
set output "data_results/pressure.png"

set key inside top left

set xlabel("Density (g/cm^3)")
set ylabel("Pressure (Pa)")

plot "data_results/enegies.dat" u 1:7 t "Total pressure" with linespoints lc rgb 'black' pt 7 ps 1,\
     "data_results/enegies.dat" u 1:6 t "Kinetic pressure" with linespoints lc rgb 'red' pt 7 ps 1,\
     "data_results/enegies.dat" u 1:5 t "Potential pressure" with linespoints lc rgb 'blue' pt 7 ps 1
     
     
