 set term postscript eps color enhanced

set xlabel 'U'
set ylabel '<M>'

set output 'caloric.eps'

pl for[ii=0:9] 'DrivenHMF_w1-10_nh20.dat' every:::ii::ii u 1:4 w lp pt 15 lw 1
