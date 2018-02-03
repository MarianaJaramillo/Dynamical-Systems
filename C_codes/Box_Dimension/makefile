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
#	-O0 : Este nivel desconecta por completo la optimización.
#	-O1 : El nivel de optimización más básico.
#	-O2 : Un paso delante de -O1.
#             Es el nivel recomendado de optimización.
#	-O3 : El nivel más alto de optimización posible.
#
#============================================================
#OPT=-O3

#============================================================
# LDFLAGS:
#	Opciones de linkado
#
# -g: activa el depurador gdb
#
# -o: Hace el linkado de los archivos .o y genera el ejecutable
#
# MATHFLAGS:
# -lm: Hace el linkado de la libreria math.h
#      (debe ir al final del comando de compilacion)
#
# GSLFLAGS=-lgsl -lgslcblas
#	Hace el linkado de las librerias de GSL
# Ej:
#
#	$ gcc OPT -Wall -g -o ejecutable.x *.o -lgsl -lgslcblas -lm
#============================================================
LDFLAGS=-Wall -g -o
MATHFLAGS=-lm
GSLFLAGS=-lgsl -lgslcblas


#============================================================
#	Box dimension
#============================================================
SOURCES_box=Henon_Box_Dimension.c Fractal_Dimension.c Dynamical_Systems_Structs.c Matrix-Vector_methods.c
OBJECTS_box=$(SOURCES_box:.c=.o)

box:  $(OBJECTS_box)
	$(CC) $(OPT) $(LDFLAGS) Henon_Box_Dimension.x $(OBJECTS_box) $(GSLFLAGS) $(MATHFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) $< -o $@


#============================================================
#	Clean ALL
#============================================================
clean:
	rm -r *.o *.x

#============================================================
#	Clean Error - Output files
#============================================================
errclean:
	rm -r *.err* *.out*
