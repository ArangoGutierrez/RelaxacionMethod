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
# -D SAVE=1 Save the boolean matrix in Maps directory as table format
# -D SAVE=2 Save the boolean matrix in Maps directory as marix format
# -D SAVE=3 Save the boolean matrix in Maps directory in both formats matrix and table
# -D ITERATIONS=<int> Cellular automata number of generations
# -D TIME to take time 
BFLAGS = -D ITERATIONS=1000 -D SAVE=3 -D TIME

runbool: clean compilebool
	@./$(TARGET_1)

compilebool:
	@$(COMPILER) $(TARGET_1)$(EXT) -o $(TARGET_1) $(CFLAGS) $(BFLAGS)

# The Boolmap should be printed as matrix format 
# Compile the program with -D SAVE=2 and run it.
# Then you can make the following command
plotbool:
	@gnuplot -p -e "filename='../Maps/boolmap_matrix.dat'" gnuplot.sh

clean:
	@rm -rf $(TARGET_1) $ $(TARGET_2) *.txt ../Maps/boolmap.dat


runrelaxation: compilerelaxation

