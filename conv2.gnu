#
#   GNUPLOT v3.6 beta multiplot script file
#
plot 'conv.dat' using 0:3 with line title  'r',\
     'conv.dat' using 0:4 with line title 'ru',\
     'conv.dat' using 0:5 with line title 'rv',\
     'conv.dat' using 0:6 with line title 'Et'
#==========================================
pause -1 "Hit return to continue"

