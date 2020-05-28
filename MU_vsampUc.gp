#code to plot amp vs max dev. for different values of Lambda for many body Kapitza pendulum

set term postscript eps color enhanced
set output 'MU_vsampUc.eps'

set multiplot
set size 1.0,1
set tic font "Times Roman, 25"
unset key 
set xtics border offset   0,-0.25,0
set ytics border offset   -0.25,0.0,0
set xrange[0.15:1]
set yrange[0:0.5]
set xlabel 'amplitude (a)' font "Times Roman,25" offset 0,-1,0
set ylabel ' U_c'  font "Times Roman,25" offset -1,0,0
set bmargin 5
set lmargin 10
set tmargin 2
set border linewidth 3
pl  'table_criticalpoint.txt'  u 1:2 w lp pt 4 lw 3 ps 1.5 lc "blue"

set origin .48, .1
set size .5,.55

#set  key box
set key font "Time Roman,13"
set key spacing 1.3
#set key box width 1
#set key box vertical width 2 height .5 maxcols 1 spacing 2

set tic font "Times Roman ,15"
set xtics border offset   0,-0.25,0
set ytics border offset   0,0.0,0
set border linewidth 1.5

set xrange[0.0:0.75]
set yrange[0:1]

set xlabel "U" font "Time Roman,20" offset   0,-.25,0
set ylabel "M" font "Time Roman,20"  offset 1.5,0.0,0


 pl  'amp0.2.dat'  u 1:4 w l lw 1 title "a = 0.2"
repl  'amp0.4.dat'  u 1:4 w l lw 1 title "a = 0.4"
repl  'amp0.5.dat'  u 1:4 w l lw 1 title "a = 0.5"
repl  'amp0.6.dat'  u 1:4 w l lw 1 title "a = 0.6"
repl  'amp0.7.dat'  u 1:4 w l lw 1 title "a = 0.7"
repl  'amp0.8.dat'  u 1:4 w l lw 1 title "a = 0.8"
repl  'amp0.9.dat'  u 1:4 w l  lw 3 title "a = 0.9"

