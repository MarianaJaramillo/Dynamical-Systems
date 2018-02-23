//
// Copyright (c) 2018 by Camilo-HG. All Rights Reserved.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// GSL libraries for handling FFT
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

// Library for handling vectors and matrices
#include "Matrix-Vector_methods.h"

// Library with defined structs
#include "Dynamical_Systems_Structs.h"

#ifndef PERIODOGRAM  // header guard. Unique identifier for each header
#define PERIODOGRAM

/*************************************************************
  Methods and Routines
*************************************************************/

void DynamicalSystem2array(Dynamical_System * sys,
                          const int coordinate,
                          double **data,
                          int *Ndata);

void Periodogram(double *data,
                const int Ndata,
                const double Delta,
                double **frecuencies,
                double **powers,
                int *Nfrecuencies);

#endif // PERIODOGRAM
