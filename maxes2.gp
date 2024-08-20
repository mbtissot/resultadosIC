# Set the terminal type and output file
# set terminal pngcairo# enhanced size 800,600

set term postscript eps color blacktext "Helvetica" 14
set output 'maxes2.eps'

# Set titles and labels
set title "Maximos: Peak2"
set xlabel "Tempo"
set ylabel "Peak2"

#set persist

# Set grid and legend
set grid
set key outside

#set logscale y
set xrange [3000:5000]
set yrange [0.0004:0.0022]

# Plot each file with a different line style
plot "0.85_maxes.txt"  u 1:3 w l, \
     "0.90_maxes.txt"  u 1:3 w l, \
     "0.95_maxes.txt" u 1:3 w l, \
     "1.00_maxes.txt"  u 1:3 w l, \
     "1.05_maxes.txt" u 1:3 w l, \
     "1.10_maxes.txt"  u 1:3 w l, \
     "1.15_maxes.txt"  u 1:3 w l



set xlabel 'u_{/Symbol \136}' font ",22";
set ylabel 'u_{/Symbol \174\174}' font ",22";


#pause mouse close
