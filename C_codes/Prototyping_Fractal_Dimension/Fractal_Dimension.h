#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

// GSL library for statistics
#include <gsl/gsl_statistics.h>

// Library for handling vectors and matrices
#include "Matrix-Vector_methods.h"

// Library with defined structs
#include "Dynamical_Systems_Structs.h"

/*************************************************************
  Methods and Routines
*************************************************************/

double logbase(double y, int b);

void Box_Dimension(const Dynamical_System * sys,
                  double *sys_Fractal_Dimension,
                  double *sys_Fractal_Dimension_SD,
                  const char *filename,
                  const double epsilon,
                  const int max_exponent,
                  const double tolerance_pct,
                  const int output_option,
                  const int print_time_option,
                  const char *time_filename);
