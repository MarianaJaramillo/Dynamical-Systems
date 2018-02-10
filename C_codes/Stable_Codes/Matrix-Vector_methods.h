//
// Copyright (c) 2018 by Camilo-HG. All Rights Reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#ifndef MATRIX_VECTOR  // header guard. Unique identifier for each header
#define MATRIX_VECTOR

/*************************************************************
  Methods and Routines
*************************************************************/

void Vector_copy(double *a,
                double *a_asst,
                int n);

void Matrix_copy(double **a,
                double **a_asst,
                int k,
                int l);

void Vector_addition(double alpha,
                    double *a,
                    double beta,
                    double *b,
                    double *c,
                    int n);

void Vector_subtraction(double alpha,
                        double *a,
                        double beta,
                        double *b,
                        double *c,
                        int n);

void Vector_Scalar_product(double *a,
                          double scalar,
                          int n,
                          double *c);

void Matrix_Scalar_product(double **a,
                          double escalar,
                          int m,
                          int n,
                          double **c);

void Matrices_product(double **a,
                      double **b,
                      int m,
                      int l,
                      int n,
                      double **c);

double Dot_product(double *a,
                  double *b,
                  int n);

double Norm(double *a,
            int n);

#endif // MATRIX_VECTOR
