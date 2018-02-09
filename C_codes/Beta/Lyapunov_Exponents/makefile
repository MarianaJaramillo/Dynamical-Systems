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
#	Exponentes de Lyapunov
#============================================================
SOURCES_lyapunov=Lyapunov_exponents.c metodos_matrices.c
OBJECTS_lyapunov=$(SOURCES_lyapunov:.c=.o)

lyapunov:  $(OBJECTS_lyapunov)
	$(CC) $(OPT) $(LDFLAGS) Lyapunov_exponents.x $(OBJECTS_lyapunov) $(GSLFLAGS) $(MATHFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) $< -o $@


#============================================================
#	Clean ALL
#============================================================
clean:
	rm -r *.o *.x
