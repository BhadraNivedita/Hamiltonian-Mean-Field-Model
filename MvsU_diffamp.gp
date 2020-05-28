#code to plot amp vs max dev. for different values of Lambda for many body Kapitza pendulum

set term postscript eps color enhanced

set xlabel "initial energy (U)" font "Time Roman,30" offset   0,-1,0
set ylabel "T" font "Time Roman,30"  offset -0.5,0.0,0

#unset key
set key box
set key font "Time Roman,20"
set key spacing 1.25
set key box width .7
set tic font "Times Roman,20"
set xtics border offset   0,-0.5,0
set ytics border offset   -0.5,0.0,0
set border linewidth 3
set key at 0.65,0.95

set xrange[0:0.75]
set yrange[0:1]

set bmargin 6
set lmargin 10
set tmargin 2
pl  'amp0.1.dat'  u 1:4 w lp pt 6 lw 3 title "a = 0.1"
repl  'amp0.2.dat'  u 1:4 w lp pt 6 lw 3 title "a = 0.2"
repl  'amp0.3.dat'  u 1:4 w lp pt 6 lw 3 title "a = 0.3"
repl  'amp0.4.dat'  u 1:4 w lp pt 6 lw 3 title "a = 0.4"
repl  'amp0.5.dat'  u 1:4 w lp pt 6 lw 3 title "a = 0.5"
repl  'amp0.6.dat'  u 1:4 w lp pt 6 lw 3 title "a = 0.6"
repl  'amp0.7.dat'  u 1:4 w lp pt 6 lw 3 title "a = 0.7"
repl  'amp0.8.dat'  u 1:4 w lp pt 6 lw 3 title "a = 0.8"
 set output 'MvsU_diffamp.eps'
repl  'amp0.9.dat'  u 1:4 w lp pt 6 lw 3 title "a = 0.9"
