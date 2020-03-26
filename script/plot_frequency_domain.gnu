inputfile = 'fftoutput.txt'

set title "Frequency Domain"
set ylabel "Amplitude"
set xlabel "Frequency [Hz]"

stats inputfile

plot \
[:STATS_records/2]inputfile using :(abs($2)/60001)
