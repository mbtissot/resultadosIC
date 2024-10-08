# Set the terminal type and output file
# set terminal pngcairo# enhanced size 800,600

# Set titles and labels
set title "Maximos: Peak1"
set xlabel "Tempo"
set ylabel "Peak1"

#set persist

# Set grid and legend
set grid
#set key outside

# Plot each file with a different line style
plot "0.85_maxes.txt"  u 1:2 w l, \
     "0.90_maxes.txt"  u 1:2 w l, \
     "0.95_maxes.txt" u 1:2 w l, \
     "1.00_maxes.txt"  u 1:2 w l, \
     "1.05_maxes.txt" u 1:2 w l, \
     "1.10_maxes.txt"  u 1:2 w l, \
     "1.15_maxes.txt"  u 1:2 w l

pause mouse close
