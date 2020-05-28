 set term postscript eps color enhanced

set xlabel 'U'
set ylabel '<M>'

unset key
set output 'MvsUamp0.9.eps'

pl  'amp0.9.dat'  u 1:4 w lp pt 15 lw 1 title "w=10"

