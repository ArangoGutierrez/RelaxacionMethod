# FILE
TARGET_1 = boolmatrix
TARGET_2 = rlxmmthd
EXT = .c
COMPILER = pgc++

# CFLAGS
# -lm Math library
# -O3 Compiler optimization
CFLAGS = -lm -O3

# BOOL MATRIX FLAGS
# -D SAVELAST=1 Save the boolean matrix in Maps directory as table format
# -D SAVELAST=2 Save the boolean matrix in Maps directory as marix format
# -D SAVELAST=3 Save the boolean matrix in Maps directory in both formats matrix and table
# -D SAVEINIT=1 Save the bool matrix contour as table format
# -D SAVEINIT=2 Save the bool matrix contour as matrix format
# -D SAVEINIT=3 Save the bool matrix contour in both formats
# -D ITERATIONS=<int> Cellular automata number of generations
# -D CONTOURS=1..3 Read a number of contour files
# -D TIME to take time 
BFLAGS = -D ITERATIONS=300 -D CONTOURS=3 -D SAVELAST=3 -D SAVEINIT=3 -D TIME

# RELAXATION METHOD FLAGS
# -D SAVEINIT Save the temperature matrix initial state
# -D SAVEALL Save all generations of the temperature matrix
# -D SAVELAST Save the temperature matrix last state
# -D GENERATIONS=<int> The system evolves until it reaches a given number of generations
# -D GENERATIONS=-1 The system evolves until it reaches a steady state (convergence)
# -D TIME to take time 
RGENERATIONS = 300
RFLAGS = -D GENERATIONS=$(RGENERATIONS) -D SAVEINIT -D SAVELAST -D TIME

# BOOL MATRIX 
runbool: clean compilebool
	@./$(TARGET_1)

compilebool:
	@$(COMPILER) $(TARGET_1)$(EXT) -o $(TARGET_1) $(CFLAGS) $(BFLAGS)

# The Boolmap should be printed as matrix format 
# Compile the program with -D SAVE=2 and run it.
# Then you can make the following command
plotbool:
	@gnuplot -p -e "filename='../OutputData/boolmap_matrix.dat'" ../OutputData/bool_gnuplot.gnu

plotboolcontour:
	@gnuplot -p -e "filename='../OutputData/boolmap_contour_matrix.dat'" ../OutputData/bool_gnuplot.gnu

#plotboolheatmap:
#	@gnuplot -p -e "filename='../OutputData/boolmap_table.dat'" ../OutputData/heatmap_gnuplot.gnu


# RELAXATION METHOD
runrelaxation: runbool compilerelaxation
	@./$(TARGET_2)


compilerelaxation:
	@$(COMPILER) $(TARGET_2)$(EXT) -o $(TARGET_2) $(CFLAGS) $(RFLAGS)

plotrlxinitheatmap:
	@gnuplot -p -e "filename='../OutputData/RlxMthd_v1.0_1.dat'" ../OutputData/heatmap_gnuplot.gnu

plotrlxlastheatmap:
	@gnuplot -p -e "filename='../OutputData/RlxMthd_v1.0_$(RGENERATIONS).dat'" ../OutputData/heatmap_gnuplot.gnu

createrelaxationimages:
	@../OutputData/heatmap_images.sh $(RGENERATIONS)

# Remove boolmatrix binary
# Remove rlxmmthd binary
# Remove boolmap_matrix and boolmap_table
# Remove time_bool and time_relaxation
# Remove generation files
clean:
	@rm -f $(TARGET_1) $(TARGET_2) ../OutputData/boolmap*.dat ../OutputData/RlxMthd_v1.0_*.dat ../OutputData/*.png
