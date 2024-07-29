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
plot "0.9_ratios.txt"  w l, \
     "0.95_ratios.txt" w l, \
     "1.0_ratios.txt"  w l, \
     "1.05_ratios.txt" w l, \
     "1.1_ratios.txt"  w l

pause mouse close