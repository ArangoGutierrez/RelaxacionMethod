# FILE
TARGET_1 = boolmatrix
TARGET_2 = rlxmmthd
EXT = .c
COMPILER = g++

# CFLAGS
# -lm Math library
# -O3 Compiler optimization
CFLAGS = -lm -O3 

# BOOL MATRIX FLAGS
# -D SAVE=1 Save the boolean matrix in Maps directory as table format
# -D SAVE=2 Save the boolean matrix in Maps directory as marix format
# -D SAVE=3 Save the boolean matrix in Maps directory in both formats matrix and table
# -D ITERATIONS=<int> Cellular automata number of generations
# -D TIME to take time 
BFLAGS = -D ITERATIONS=1000 -D SAVE=3 -D TIME

# RELAXATION METHOD FLAGS
# -D SAVEINIT Save the temperature matrix initial state
# -D SAVEALL Save all generations of the temperature matrix
# -D SAVELAST Save the temperature matrix last state
# -D GENERATIONS=<int> The system evolves until it reaches a given number of generations
# -D GENERATIONS=0 The system evolves until it reaches a steady state (convergence)
# -D TIME to take time 
RGENERATIONS = 1
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
	@gnuplot -p -e "filename='../Maps/boolmap_matrix.dat'" bool_gnuplot.sh


# RELAXATION METHOD
runrelaxation: runbool compilerelaxation
	@./$(TARGET_2)

compilerelaxation:
	@$(COMPILER) $(TARGET_2)$(EXT) -o $(TARGET_2) $(CFLAGS) $(RFLAGS)

plotheatmap:
	@./heatmap_gnuplot.sh $(RGENERATIONS)

# Remove boolmatrix binary
# Remove rlxmmthd binary
# Remove boolmap_matrix and boolmap_table
# Remove time_bool and time_relaxation
# Remove generation files
clean:
	@rm -f $(TARGET_1) $(TARGET_2) ../Maps/boolmap*.dat time*.txt ../OutputData/RlxMthd_v1.0_*.dat *.png
