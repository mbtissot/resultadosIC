# Set the terminal type and output file
# set terminal pngcairo# enhanced size 800,600

set term postscript eps color blacktext "Helvetica" 14
set output 'ratioslog.eps'

# Set titles and labels
set title "Ratios: Peak2/Peak1"
set xlabel "Tempo"
set ylabel "P_2/P_1"

#set persist

set logscale y
set xrange [3000:5000]
set yrange [0.025:0.08]

# Set grid and legend
set grid
set key outside

# Plot each file with a different line style
plot "0.85_ratios.txt" w l, \
     "0.90_ratios.txt"  w l, \
     "0.95_ratios.txt" w l, \
     "1.00_ratios.txt"  w l, \
     "1.05_ratios.txt" w l, \
     "1.10_ratios.txt"  w l, \
     "1.15_ratios.txt"  w l

#pause mouse close
