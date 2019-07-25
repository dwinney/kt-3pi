set term pdf color enhanced font "Helvetica"
set auto
set size square
set border 15 lw 1

set xlabel "√s [GeV]"
set xrange [0.3:1.]

set xtics 0.2

set out 'std_once_vs_poly_1.pdf'
set multiplot
set title ""
set key left
p "./std_one_subtraction.dat" u 1:4 t "" with l lw 2 lc "dark-grey", \
  "./fit_results_1.dat" u 1:2 t " " with l dt 2 lw 2 lc "black", \
  "./fit_results_2.dat" u 1:2 t " " with l dt 3 lw 2 lc "black", \
  "./fit_results_3.dat" u 1:2 t " " with l dt 4 lw 2 lc "black", \
