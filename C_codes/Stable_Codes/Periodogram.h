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

//============================================================

void Periodogram(const double data[],
                const int Ndata,
                const double Delta,
                double *frecuencies,
                double *powers,
                int *Nfrecuencies);
