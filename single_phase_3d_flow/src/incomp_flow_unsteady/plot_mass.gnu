set terminal pngcairo size 600,400 enhanced font 'Verdana,12' 
set for [i=1:9] linetype i dashtype i
set xlabel "Time (s)"
set ylabel "Mass (kg)" offset 1,0
set output 'mass.png'
set mxtics 5 
set mytics 5 
#set tic scale 1.2
p'time_massflow.txt' u 2:($3) w l lw 2 lt 1 lc rgb "blue" title "Total mass in",\
 'time_massflow.txt' u 2:($4) w l lw 2 lt 1 lc rgb "red" title "total mass out"
