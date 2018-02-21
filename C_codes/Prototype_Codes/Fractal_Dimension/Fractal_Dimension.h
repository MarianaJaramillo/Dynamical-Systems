//
// Copyright (c) 2018 by Camilo-HG. All Rights Reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

// GSL library for statistics
#include <gsl/gsl_statistics.h>

// GSL library for sorting vectors
#include <gsl/gsl_sort_double.h>

// Library for handling vectors and matrices
#include "Matrix-Vector_methods.h"

// Library with defined structs
#include "Dynamical_Systems_Structs.h"

#ifndef FRACTAL_DIMENSION  // header guard. Unique identifier for each header
#define FRACTAL_DIMENSION

/*************************************************************
  Define structs
*************************************************************/
typedef struct Hyperbox {
  /*============================================================
    STRUCT: Hyperbox

      box[sys.dimension][2]
      ==========================


      points[Npoints][sys.dimension]
      ==========================

        **points: Stores system's points belonging

        Npoints: Gives us information about how
                many points has the time series.

        dimension: Gives us information about the
                  dimensionality of the system.

        Thus:

          points[i][j] ---> corresponds to the j-th component
                            of the i-th point.
  ============================================================*/
  double **box;
  double **point;
  size_t Npoints;
} Hyperbox;


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

#endif // FRACTAL_DIMENSION
