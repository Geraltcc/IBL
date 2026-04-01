set terminal pngcairo size 1200,800 enhanced font 'Arial,14'
set output 'data/ue_grad.png'

set xlabel '{x}'
set ylabel '{ue_{grad}}'
set title 'ue_{grad}'
set grid
set key top right

plot 'data/ue_grad.dat' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'present', \
 