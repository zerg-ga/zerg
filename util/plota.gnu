
set style line 1 lt 1 lc rgb "yellow" lw 2
set style line 2 lt 1 lc rgb "purple" lw 5
set style line 3 lt 1 lc rgb "red" lw 6
set style line 4 lt 1 lc rgb "black" lw 7
set style line 5 lt 1 lc rgb "blue" lw 6
set style line 6 lt 1 lc rgb "cyan" lw 3 
set style line 7 lt 1 lc rgb "green" lw 2

plot "meanRun.txt" using 1:2 smooth acsplines title 'Immigration' ls 1, \
     "meanRun.txt" using 1:3 smooth acsplines title 'Sphere Crossover' ls 2, \
     "meanRun.txt" using 1:4 smooth acsplines title 'Geometric Center Displacement' ls 3, \
     "meanRun.txt" using 1:5 smooth acsplines title 'TwistOperator' ls 4, \
     "meanRun.txt" using 1:6 smooth acsplines title 'Cut and Splice' ls 5, \
     "meanRun.txt" using 1:7 smooth acsplines title 'Angular Operator' ls 6, \
     "meanRun.txt" using 1:8 smooth acsplines title 'Surface Angular Operator' ls 7



