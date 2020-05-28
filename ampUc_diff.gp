#code to plot amp vs max dev. for different values of Lambda for many body Kapitza pendulum

set term postscript eps color enhanced

set xlabel "amplitude ( a )" font "Time Roman,30" offset   0,-1,0
set ylabel "Critical point ( U_c )" font "Time Roman,30"  offset -0.5,0.0,0

#unset key
set key box
set key font "Time Roman,25"
set key spacing 1.25
set key box width .7
set key at 0.3,.65
set tic font "Times Roman,25"
set xtics border offset   0,-0.5,0
set ytics border offset   -0.5,0.0,0
set border linewidth 3

set xrange[0:1]
set yrange[0:0.75]

set bmargin 6
set lmargin 10
set tmargin 2


 pl for[ii=0:0] 'table_criticalpoint_w0.txt'  u 1:2 w lp pt 4  lw 3 lc "black" title "w=0"
  set output 'ampUc_diff.eps'
  repl for[ii=0:0] 'table_criticalpoint.txt'  u 1:2 w l dt 4 ps 1.5 lw 3 lc "blue" title "w=10"
