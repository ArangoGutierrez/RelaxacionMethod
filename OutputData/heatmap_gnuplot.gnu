reset
set view map
set pm3d interpolate 0,0
set dgrid3d

#set cbrange [0.0000:30.0000]
#set palette rgb 7,5,15
#splot filename using 1:2:3 w image

#set cbrange [0:30]
#set palette defined (0 "blue", 30 "red")
#splot filename using 1:2:3 with pm3d

#set cbrange [0:6]
set palette rgb 7,5,15
splot filename using 1:2:3 with pm3d
