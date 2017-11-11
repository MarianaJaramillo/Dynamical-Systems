#============================================================
#		MY MAKEFILE
#============================================================

#============================================================
# Compilador a usar
#============================================================
CC=gcc

#============================================================
# CFLAGS:
# 	Opciones de compilacion
#
# -Wall: Muestra todas las advertencias
#
# -g: Activa el depurador gbd
#
# -c: Compila y genera los archivos .o
#
# Ej:
#
#	$ gcc -Wall -g -c *.c -lm
#============================================================
CFLAGS=-Wall -g -c

#============================================================
# OPT:
#	Opciones de optimizacion
#
# -O: Variable que controla el nivel de
#     optimizacion de todo el codigo
#     hecho por el propio compilador.
#
#     Niveles de optimizacion:
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
#	Periodograma
#============================================================
SOURCES_periodograma=Periodograma.c
OBJECTS_periodograma=$(SOURCES_periodograma:.c=.o)

periodograma:  $(OBJECTS_periodograma)
	$(CC) $(OPT) $(LDFLAGS) Periodograma.x $(OBJECTS_periodograma) $(GSLFLAGS) $(MATHFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) $< -o $@


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

#============================================================
#	Clean Error - Output files
#============================================================
errclean:
	rm -r *.err* *.out*
