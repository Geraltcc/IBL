set terminal pngcairo size 1200,800 enhanced font 'Arial,14'
set output 'data/ue_vec.png'

set xlabel '{/Symbol x} (arc length from stagnation)'
set ylabel '{/Symbol ue} (surface fractions)'
set title 'IBL Surface Fraction Coeff'
set grid
set key top right
set xrange[0:1]

plot 'data/ue_vec.dat' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'present'