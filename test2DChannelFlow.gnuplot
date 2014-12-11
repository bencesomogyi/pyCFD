set datafile separator ","
set title "Residuals"
set logscale y
set xrange[:]
set yrange [1e-4:]
set xlabel "Time [s]"
set ylabel "Residual"
set grid xtics
set grid ytics
plot "_OUTPUT/residuals_p.csv" using 1:2 title "res_p_corr" with lines,\
     "_OUTPUT/residuals_U.csv" using 1:2 title "res_U_X" with lines,\
     "_OUTPUT/residuals_U.csv" using 1:3 title "res_U_Y" with lines
	 
pause 1
reread
