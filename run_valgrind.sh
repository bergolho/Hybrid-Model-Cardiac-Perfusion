#!/bin/bash
PNAME="./gadolino"
ARGS="0.001 100 meshes/arvore_1.msh meshes/arvore_2.msh meshes/arvore_3.msh 100 meshes/vent100x100.dat 8"

if [ ! -f $PNAME ]; then
	./recompile_project.sh
fi

valgrind --leak-check=full --show-leak-kinds=all $PNAME $ARGS
