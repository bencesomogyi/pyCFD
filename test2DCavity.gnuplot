set datafile separator ","
set term x11 0
set title "Residuals"
set logscale y
set xlabel "Time [s]"
set ylabel "Residual"
set grid xtics
set grid ytics
plot "_OUTPUT/residuals_p.csv" using 1:2 title "res_p_corr" with lines,\
	 "_OUTPUT/residuals_U.csv" using 1:2 title "res_U_X" with lines,\
	 "_OUTPUT/residuals_U.csv" using 1:3 title "res_U_Y" with lines
	 
set term x11 1
set title "Monitor p"
unset logscale
set ylabel "p"
plot "_OUTPUT/residuals_p.csv" using 1:3 title "p" with lines
	 
pause 1
reread
