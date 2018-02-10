//
// Copyright (c) 2018 by Camilo-HG. All Rights Reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>


/*************************************************************
 METODOS Y RUTINAS
*************************************************************/

void Copiar_Vectores(double *a,
                    double *a_aux,
                    int n);

void Copiar_Matrices(double **a,
                    double **a_aux,
                    int k,
                    int l);

void Sumar_Vectores(double alpha,
                    double *a,
                    double beta,
                    double *b,
                    double *c,
                    int n);

void Restar_Vectores(double alpha,
                    double *a,
                    double beta,
                    double *b,
                    double *c,
                    int n);

void Producto_Vector_Escalar(double *a,
                            double escalar,
                            int n,
                            double *c);

void Producto_Matriz_Escalar(double **a,
                            double escalar,
                            int m,
                            int n,
                            double **c);

void Producto_Matrices(double **a,
                      double **b,
                      int m,
                      int l,
                      int n,
                      double **c);

double Producto_Interno(double *a,
                      double *b,
                      int n);

double Norma(double *a,
            int n);
