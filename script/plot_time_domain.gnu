set terminal png
set title "Time Domain"
set ylabel "Amplitude"
set xlabel "Time [s]"
set output 'time_domain.png'

plot \
"dumpdata.txt"
