# Set the terminal type and output file
# set terminal pngcairo# enhanced size 800,600

# Set titles and labels
set title "Maximos: Peak2"
set xlabel "Tempo"
set ylabel "Peak12"

#set persist

# Set grid and legend
set grid
#set key outside

# Plot each file with a different line style
plot "0.9_maxes.txt"  u 1:3 w l, \
     "0.95_maxes.txt" u 1:3 w l, \
     "1.0_maxes.txt"  u 1:3 w l, \
     "1.05_maxes.txt" u 1:3 w l, \
     "1.1_maxes.txt"  u 1:3 w l

pause mouse close