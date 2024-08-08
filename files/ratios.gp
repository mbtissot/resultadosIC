# Set the terminal type and output file
# set terminal pngcairo# enhanced size 800,600

# Set titles and labels
set title "Ratios: Peak2/Peak1"
set xlabel "Tempo"
set ylabel "P_2/P_1"

#set persist

# Set grid and legend
set grid
#set key outside

# Plot each file with a different line style
plot "0.80_ratios.txt" u 1:2 w l, \
     "0.85_ratios.txt" u 1:2 w l, \
     "0.90_ratios.txt" u 1:2 w l, \
     "0.95_ratios.txt" u 1:2 w l, \
     "1.00_ratios.txt" u 1:2 w l, \
     "1.05_ratios.txt" u 1:2 w l, \
     "1.10_ratios.txt" u 1:2 w l,\
     "1.15_ratios.txt" u 1:2 w l

pause mouse close