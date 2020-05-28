#code to show fluctuation in integration value of abs(cos(wt))


set term postscript eps color enhanced
set output 'datafluctuation.eps'


set multiplot
set size 1.0,1

unset key
set xtics font ",30"  
set ytics font ",30"  
set xtics border offset   0,-0.5,0
set ytics border offset   -0.5,0.0,0
set border linewidth 3
set bmargin 6
set lmargin 16
#set object ellipse center .13, 0 size .4, 4
#set arrow from .1, 2.1 to screen .22, .4 front lt 3

set yrange[3.95:4.065]
set xrange[5:20]
set xlabel 'n' font ",40" offset 0,-1,0
set ylabel 'f(t)'  font ",40" offset -5,0,0
pl 'datafluctuation.dat' u 1:2 w lp lt -1

set origin .6, .68
set size .3,.28
clear
unset key
unset xrange
set yrange [0:1]
set border linewidth 2.0
set tics font ",15"
#unset grid
#unset object
#unset arrow
set xlabel '{/Symbol f}'  offset 0,0,0
set ylabel '{/Symbol e}(t)' offset 0,0,0
set bmargin 1
set tmargin 1
set lmargin 3
set rmargin 1
plot [0:2*pi]  abs(sin(x)) w lp lw 4

unset multiplot