set terminal pngcairo size 1200,800 enhanced font 'Arial,14'
set output 'data/delta.png'

set xlabel '{/Symbol x} (arc length from stagnation)'
set ylabel '{\delta^*}'
set title 'IBL'
set grid
set key top right
set xrange[0:1]
set yrange[0:0.006]

plot 'data/delta.dat' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'present', \
     'data/BL2.dat' every ::1 using 2:5 with lines lw 2 lc rgb 'red' title 'xfoil'