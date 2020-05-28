#code to plot trapping for N=10000 and w=1

set terminal postscript color eps enhanced

set key box 
set key font "Times-Roman,25"
set key spacing 1.25
set key box width 1.80

set border linewidth 3
set tics font "Times-Roman,25" 
set xtics border offset   0,-0.5,0
set ytics border offset   -0.5,0.0,0
set lmargin 12 
set bmargin 5
set xlabel "U" font "Times-Roman,30" offset 0,-1,0
set ylabel "P(U)" font  "Times-Roman,30"  offset -2,0,0
set yrange[0:1.2]
set xrange[0:1]
	#pl 'trapping_w1.dat'  w lp pt 6 linewidth 4 title "w=1"
	pl 'trapping_w0.dat'  w lp pt 4 ps 1.6 lw 4 lc "black" title "w=0"
	repl 'trapping_w10.dat'  w lp pt 3 ps 2 lw 4 lc "black"  title "w=10"
set output 'driventrapping.eps' 
	
	repl 'trapping_data20.dat'  w lp pt 6 ps 2 dt 4  lw 4 lc "black" title "w=20"
