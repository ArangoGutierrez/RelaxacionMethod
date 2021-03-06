#bin/bash
for i in $(seq 1 $1)
do
    /usr/bin/gnuplot <<__EOF
set view map
#set pm3d
# More fine picture
set pm3d interpolate 0,0
set dgrid3d
# Set Palette Heatmap 
set palette rgb 7,5,15 
set title "Relaxation Method"
set term png
set output "../OutputData/image-${i}.png"
splot "../OutputData/RlxMthd_v1.0_${i}.dat" using 1:2:3 with pm3d
pause -1
__EOF
done
