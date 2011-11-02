#
# Plot TINKER Monthly/Cummulative Downloads/Webhits
#
set output "plot.ps"
set terminal postscript color
set style data lines
set nokey
set ytics nomirror
set y2tics nomirror
set xrange [0:180]
set yrange [0:3000]
set y2range [0:80]
plot "stats" using 1:4 axes x1y1, "stats" using 1:5 axes x1y2
set terminal postscript color
set style data lines
set nokey
set ytics nomirror
set y2tics nomirror
set xrange [0:180]
set yrange [0:7000]
set y2range [0:500]
plot "stats" using 1:6 axes x1y1, "stats" using 1:7 axes x1y2
quit
