#
# Copyright (c) 2018 by Camilo-HG. All Rights Reserved.
#

#============================================================
#	MAKEFILE
#============================================================

#============================================================
# COMPILER TO USE
#============================================================
CC=gcc

#============================================================
# CFLAGS:
#		Compilation options
#
# -Wall: shows every warning
#
# -g: Activates gbd debugger
#
# -c: Compiles and generates .o files
#
# Example:
#
#	$ gcc -Wall -g -c *.c -lm
#============================================================
CFLAGS=-Wall -g -c

#============================================================
# OPT:
#		Optimization options
#
# -O: Variable which controls optimization
#     level of the whole code done by the
#     compiler itself.
#
#     Optimization levels:
#	-O0 : optimization disabled
#	-O1 : Most basic optimization level
#	-O2 : One step forward -O1.
#             This is the recommended optimization level.
#	-O3 : Higher possible optimization level.
#
#============================================================
#OPT=-O3

#============================================================
# LDFLAGS:
#	Linkage's flags
#
# -g: Activates debugging with GDB
#
# -o: Links .o objects and generates executable program
#
# MATHFLAGS:
# -lm: links to math.h
#      (it must be located at the final place in compilation statement)
#
# GSLFLAGS=-lgsl -lgslcblas
#	Links to GSL libraries
# e.g.:
#
#	$ gcc OPT -Wall -g -o executable.x *.o -lgsl -lgslcblas -lm
#============================================================
LDFLAGS=-Wall -g -o
MATHFLAGS=-lm
GSLFLAGS=-lgsl -lgslcblas


#============================================================
#	N dimensional FOR cycle
#============================================================
SOURCES=Ndim_For_Cycle.c
OBJECTS=$(SOURCES:.c=.o)

main:  $(OBJECTS)
	$(CC) $(OPT) $(LDFLAGS) Ndim_For_Cycle.x $(OBJECTS) $(GSLFLAGS) $(MATHFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) $< -o $@


#============================================================
#	Clean ALL
#============================================================
clean:
	rm -r *.o *.x
