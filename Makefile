INCLUDE_DIR = include 
SOURCE_DIR = src
OBJECT_DIR = obj
CC = g++
PROGRAM_FLAGS =  -DVTK -DOUTPUT -DPRESSURE -DRADIUS 
CFLAGS = -Wall -Wno-write-strings -std=c++11 -fopenmp -O2 -I$(INCLUDE_DIR) $(PROGRAM_FLAGS)
CC_SOURCES = $(wildcard src/*)
STRING_OBJ_AUX = $(CC_SOURCES:.cpp=.o)
STRING_OBJ = $(subst src/,,$(STRING_OBJ_AUX))
CC_OBJ = $(patsubst %,$(OBJECT_DIR)/%,$(STRING_OBJ))
PROGRAM_NAME = gadolino

all: $(PROGRAM_NAME)

$(PROGRAM_NAME): $(CC_OBJ)
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	$(CC) $(CFLAGS) -o $@ -c $< -lm

clean:
	rm -f $(OBJECT_DIR)/*.o $(PROGRAM_NAME)
	rm -f VTK/*.vtk
	rm -f *.vtk
	rm -f sol/tissue/*.vtk
	rm -f sol/gauss.txt
	rm -f sol/graph1/*.vtk
	rm -f sol/graph2/*.vtk
	rm -f sol/graph3/*.vtk
	rm -f Output/*.dat Output/*.txt Output/*.pdf

clcResults:
	rm -f VTK/*.vtk
	rm -f *.vtk
	rm -f Output/*.dat Output/*.txt Output/*.pdf

remade:
	$(MAKE) clean
	$(MAKE)

plot:
	cd Output; python makePlot.py

print-%  : ; @echo $* = $($*)
