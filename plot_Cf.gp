set terminal pngcairo size 1200,800 enhanced font 'Arial,14'
set output 'data/Cf2_plot.png'

set xlabel '{/Symbol x} (arc length from stagnation)'
set ylabel '{/Symbol Cf} (surface fractions)'
set title 'IBL Surface Fraction Coeff'
set grid
set key top left
set xrange [0:1]

plot 'data/Cf2.dat' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'present', \
     'data/BL2.dat' every ::1 using 2:7 with lines lw 2 lc rgb 'red' title 'xfoil'