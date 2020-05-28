set term postscript eps color enhanced
set output 'Tdeviation.eps'


set multiplot
set size 1.0,1

set key at 1.44,.85
set key font "Times Roman,25"
set key spacing 1.25
set key box width 1.80
set xtics font ",25"  
set ytics font ",25"  
set xtics border offset   0,-0.6,0
set ytics border offset   -0.5,0.0,0
set border linewidth 3
set bmargin 6
set lmargin 12
set yrange[0:2.5]
set xrange[0:1.5]
set xlabel 'U' font "Times Roman,25" offset 0,-2,0
set ylabel 'T'  font "Times Roman,25" offset -2,0,0
pl 'undriven_caloric.dat' u 1:2 w l dt 3  linewidth 4 lc "black" title "w=0", 'DrivenHMF_N1000_w10.dat' u 1:3 w l dt 6 lw 4 lc "blue" title "w=10",'DrivenHMF_N1000_w20.dat' u 1:3 w lp  dt 3 pt 3   lw 4 lc "red" title "w=20"


set origin .28, .6
set size .35,.35
clear
unset key
set xrange[-0.01:1.5]
set yrange [0:1.5]
set border linewidth 2.0
set tics font ",15"
set xlabel 'U'  offset 0,-0.5,0 font "Times Roman,20"
set ylabel 'deviation in T' offset 0,0,0 font "Times Roman,20"
set bmargin 1
set tmargin 1
set lmargin 3
set rmargin 1
plot 'datadeviation_T.dat' u 1:2 w lp pt 15 lw 3 lc "black"

unset multiplot
