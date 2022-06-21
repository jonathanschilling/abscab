#!/usr/bin/gnuplot

# Encoding (the one and only)
set encoding utf8

# Use the postscript eps version because it uses the chosen font for everything but symbols
# so Greek letters look really nice!
set terminal postscript eps enhanced color size 9cm,6cm font "TeXGyreHeros-Regular,20"\
	fontfile "/usr/share/texmf-dist/fonts/type1/public/tex-gyre/qhvr.pfb"

# Design configurations
set xlabel "number of polygon corners"
set ylabel "relative deviation from reference" offset 0,-0.5

set xrange [8:1.3e9]
set yrange [-16.5:-0.5]

set format y "%g"

set logscale x

set xtics out nomirror ("10^1" 1e1, \
                        "10^2" 1e2, \
                        "10^3" 1e3, \
                        "10^4" 1e4, \
                        "10^5" 1e5, \
                        "10^6" 1e6, \
                        "10^7" 1e7, \
                        "10^8" 1e8, \
                        "10^9" 1e9  )

set ytics out nomirror ("10^{-2}" -2, \
                        "10^{-4}" -4, \
                        "10^{-6}" -6, \
                        "10^{-8}" -8, \
                        "10^{-10}" -10, \
                        "10^{-12}" -12, \
                        "10^{-14}" -14, \
                        "10^{-16}" -16  )

set grid

# first only standard summation
set output "../article/img/McGreivy_convergence_1.eps"

plot 'convergenceMcGreivy_StandardSummation.dat'    u 1:2 w lp lc rgb "orange"          lw 2 title "on-loop  (+=)", \
     'convergenceMcGreivy_StandardSummation.dat'    u 1:3 w lp lc rgb "light-blue"      lw 2 title "McGreivy (+=)"

# now including Kahan-Babushka summation
set output "../article/img/McGreivy_convergence_2.eps"

plot 'convergenceMcGreivy_StandardSummation.dat'    u 1:2 w lp lc rgb "orange"          lw 2 title "on-loop  (+=)", \
     'convergenceMcGreivy_StandardSummation.dat'    u 1:3 w lp lc rgb "light-blue"      lw 2 title "McGreivy (+=)", \
     'convergenceMcGreivy_CompensatedSummation.dat' u 1:2 w lp lc rgb "red"        dt 2 lw 2 title "on-loop  (K-B)", \
     'convergenceMcGreivy_CompensatedSummation.dat' u 1:3 w lp lc rgb "blue"       dt 2 lw 2 title "McGreivy (K-B)"
