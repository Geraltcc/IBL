set terminal pngcairo size 1200,800 enhanced font 'Arial,14'
set output 'data/Re_th_plot.png'

set xlabel '{/Symbol x} (arc length from stagnation)'
set ylabel '{/Symbol Q}'
set title 'IBL'
set grid
set key top left
set xrange[0:1]

plot 'data/Re_th.dat' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'present'