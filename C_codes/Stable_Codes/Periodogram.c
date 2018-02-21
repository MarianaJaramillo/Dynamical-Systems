//
// Copyright (c) 2018 by Camilo-HG. All Rights Reserved.
//

// Library for handling vectors and matrices
#include "Periodogram.h"

//============================================================

void Periodogram(const double data[],
                const int Nx,
                const double Delta
                double *frecuencies,
                double *powers,
                const char *filename) {
  /*============================================================
    The program uses the fast Fourier transform of GSL to
    calculate the periodogram of the time series.

    - data[]: data's vector

    - *frecuencies = NULL

    - *powers = NULL

    - Nx : number of sampling points.

    - Delta : Sampling interval time

    Output:

      periodogram.dat ---> file with the frecuencies vector
                            and its corresponding powers.

                          Output format:

                          column 1: frecuencies
                          column 2: powers
  ============================================================*/

  // command line input variables
  int Nx;
  double Delta;

  // Program's internal variables
  int i;
  double *data_aux=NULL;
  double div_power, div_frecuencies;

  // Output files
  FILE *file=NULL;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /*=========================================================
    Asign a memory block to the vectors:

      double *data_aux;
      double *frecuencies;
      double *powers;

    Nx/2 is an integer division.
  =========================================================*/
  data_aux = (double *) malloc((size_t) Nx * sizeof(double));
  frecuencies = (double *) malloc((size_t) (Nx/2 + 1) * sizeof(double));
  powers  = (double *) malloc((size_t) (Nx/2 + 1) * sizeof(double));

  /*=========================================================
    Copy data to data_aux to avoid to modify data
  =========================================================*/
  Vector_copy(data,
              data_aux,
              Nx);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /*============================================================
    FFT variables
  ============================================================*/
  gsl_fft_real_wavetable *real;
  gsl_fft_real_workspace *work;

  // memory allocation
  work = gsl_fft_real_workspace_alloc(Nx);
  real = gsl_fft_real_wavetable_alloc(Nx);

  // FFT
  // Returns transformed data
  gsl_fft_real_transform(data_aux, 1, Nx, real, work);

  // Free memory
  gsl_fft_real_wavetable_free(real);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /*============================================================
    Periodogram
  ============================================================*/

  div_frecuencies = 1.0 / ( (double) Nx * Delta);
  div_power = 1.0/ (double) (Nx * Nx);

  // For frecuency 0
  frecuencies[0] = 0.0;

  powers[0] = div_power * (data_aux[0]*data_aux[0]);

  // For internal frecuencies
  for(i = 1; i < (Nx/2); i++){
    frecuencies[i] = (double) i * div_frecuencies;

    powers[i] = 2.0 * div_power
                  * (data_aux[2*i - 1]*data_aux[2*i - 1] + data_aux[2*i]*data_aux[2*i]);
  }

  // For Nyquist frecuency
  i = Nx/2;
  frecuencies[i] = 1.0 / (2.0 * Delta);

  if( Nx % 2 == 0) {
    powers[i] = div_power * (data_aux[Nx - 1] * data_aux[Nx - 1]);
  }
  else {
    powers[i] = div_power
                  * (data_aux[Nx - 2]*data_aux[Nx - 2] + data_aux[Nx - 1]*data_aux[Nx - 1]);
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /*============================================================
    Print
  ============================================================*/
  file = fopen(filename,"w");

  fprintf(file,
    "frecuencies\tpowers\n"
  );

  for(i = 0; i < (Nx/2 + 1); i++){
      fprintf(file,
        "%.16g\t%.16g\n"
        , frecuencies[i], powers[i]
      );
  }

  // Close file
  fclose(file);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // Free memory
  free(data_aux);

  //free(frecuencies);
  //free(powers);
}
