# test script, DO NOT USE IT!

inputfile = 'fftoutput.txt'

set title "Frequency Domain"
set ylabel "Amplitude"
set xlabel "Frequency [Hz]"

stats inputfile

set xrange[0:50]

plot \
inputfile using ($1):(abs($2)/STATS_records)
#[:STATS_records/2]inputfile using ($1):(abs($2)/60001)
