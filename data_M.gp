set term postscript eps color enhanced
set output 'caloricdeviation.eps'


set multiplot
set size 1.0,1

set key at .3,0.25
set key font ",20"
set key spacing 1.25
set key box width 1.80
set xtics font "Times Roman,30"  
set ytics font "Times Roman,30"  
set xtics border offset   0,-0.5,0
set ytics border offset   -0.5,0.0,0
set border linewidth 3
set bmargin 6
set lmargin 16


set yrange[0:1]
set xrange[0:1]
set xlabel 'U' font "Times Roman,25" offset 0,-2,0
set ylabel 'M'  font "Times Roman,25" offset -4,0,0
pl 'hmf1000.dat' u 1:3 w l dt 3  linewidth 4 lc "black" title "w=0" , 'nh50new.dat' u 1:4 w lp pt 6 lw 4 lc "blue" title "w=10",'w20.dat' u 1:4 w lp   dt 6   lw 4 lc "red" title "w=20




set origin .65, .65
set size .31,.31
clear
unset key
set xrange[-0.01:1.5]
set yrange [0:0.8]
set border linewidth 2.0
set tics font "Times Roman,15"
set xlabel 'U'  offset 0,-0.5,0 font "Times Roman,20"
set ylabel 'deviation in M' offset .25,0,0 font "times Roman,20"
set bmargin 1
set tmargin 2
set lmargin 1
set rmargin 1
plot 'datadeviation.dat' u 1:2 w lp pt 15 lw 3 lc "black"

unset multiplot
