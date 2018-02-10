//
// Copyright (c) 2018 by Camilo-HG. All Rights Reserved.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// Librería para manejo de matrices
#include "metodos_matrices.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  DEFINICIÓN DE FUNCIONES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*============================================================
  Description of the system.

  y = f(t,yo)

  The system is specified by giving the right-hand-side
  of the equation and possibly a jacobian function.

  Some methods require the jacobian function, which calculates
  the matrix dfdy and the vector dfdt.

  User-supplied parameter data is also present.
*/
typedef struct Dynamical_System {
void (*function) (size_t t, const double yo[], double f[], double *params);
size_t dimension;
double *yo;
double *params;
} Dynamical_System;


//============================================================

/* Function evaluation macros */
#define FN_EVAL(S, t, yo, f)  (*((S)->function))(t,yo,f,(S)->params)


//============================================================

void System_Map (size_t t,
                const double yo[],
                double f[],
                double *params) {
  /*============================================================
    Retorna por referencia la solución al mapa de Hénon

      Hénon(x,y) = f(x,y) = (a - x^2 + b y , x)

      Hénon map doesn't depend on t.
  ============================================================*/

(void)(t); /* avoid unused parameter warnding */
double a = *(params+0);
double b = *(params+1);

f[0] = a - yo[0]*yo[0] + b * yo[1];
f[1] = yo[0];
}


//============================================================

void Largest_Lyapunov_Exponent(const Dynamical_System * sys,
                              const double DRo,
                              const size_t transient,
                              const size_t Nmax,
                              double *Largest_Lyapunov) {
  int i;
  double t = 0.0, Df_norm;
  double *R=NULL, *R1=NULL, *f=NULL, *f1=NULL, *Df=NULL;

  *Largest_Lyapunov = 0.0;

  R = (double *) malloc((size_t) sys->dimension * sizeof(double));
  R1 = (double *) malloc((size_t) sys->dimension * sizeof(double));
  f = (double *) malloc((size_t) sys->dimension * sizeof(double));
  f1 = (double *) malloc((size_t) sys->dimension * sizeof(double));
  Df = (double *) malloc((size_t) sys->dimension * sizeof(double));

  Copiar_Vectores(sys->yo, R, sys->dimension);

  for(i = 0; i < transient; i++) {
    FN_EVAL(sys, t, R, f);
    Copiar_Vectores(f, R, sys->dimension);
  }

  Producto_Vector_Escalar(R,
                          1.0 + DRo,
                          sys->dimension,
                          R1);

  for(i=1; i < (Nmax +1); i++) {
    FN_EVAL(sys, t, R, f);
    FN_EVAL(sys, t, R1, f1);

    Restar_Vectores(1.0,
                    f1,
                    1.0,
                    f,
                    Df,
                    sys->dimension);

    Df_norm = Norma(Df,
                    sys->dimension);

    *Largest_Lyapunov += log(Df_norm/DRo);

    Copiar_Vectores(f, R, sys->dimension);

    Sumar_Vectores(1.0,
                  f,
                  DRo/Df_norm,
                  Df,
                  R1,
                  sys->dimension);
  }
  *Largest_Lyapunov = *Largest_Lyapunov/Nmax;
}

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

int main(int argc, char *argv[]) {
  /*============================================================
    Largest Lyapunov exponent
  ============================================================

    El programa calcula el exponente más grande de Lyapunov
    para el mapa de Henón.

    El programa genera el siguiente archivo de salida:

      largest_lyapunov.dat ---> archivo con el exponente más largo
                                de Lyapunov calculado hasta
                                Nmax iteraciones.

    EJECUCION DEL PROGRAMA

      ./Largest_Lyapunov_exponent.x x0 y0 a b Nmax DRo

    Donde:

      - x0 : condición inicial de x

      - y0 : condición inicial de y

      - a : valor de a

      - b : valor de b

      - Nmax : Número máximo de iteraciones

      - DRo : Separación de la segunda condición inicial
  ============================================================*/

  // Variables
  int Nmax;
  double DRo, transient;

  // Variables internas del programa
  int dim=2;
  double yo[2], params[2];
  double Lyapunov_exponent;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  yo[0] = 0.5;
  yo[1] = 0.5;
  params[0] = 1.4;
  params[1] = 0.3;
  Nmax = 1000;
  DRo  = 1e-10;
  transient = 1000;

  Dynamical_System sys = {System_Map, dim, yo, params};

  Largest_Lyapunov_Exponent(&sys,
                            DRo,
                            transient,
                            Nmax,
                            &Lyapunov_exponent);

  printf("%.16lf\n", Lyapunov_exponent);

  return 0;
}
