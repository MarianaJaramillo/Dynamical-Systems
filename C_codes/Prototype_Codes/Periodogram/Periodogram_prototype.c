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

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FUNCTION DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

//============================================================

double function(double x) {
  double fx;

  fx = sin(x);

  return fx;
}

//============================================================

void Periodogram(double data[],
                  double frecuencies[],
                  double powers[],
                  int Nx,
                  double Delta) {
  int i;

  double div_power, div_frecuencies;

  div_frecuencies = 1.0 / ( (double) Nx * Delta);
  div_power = 1.0/ (double) (Nx * Nx);

  // For frecuency 0
  frecuencies[0] = 0.0;

  powers[0] = div_power * (data[0]*data[0]);

  // For internal frecuencies
  for(i = 1; i < (Nx/2); i++){
    frecuencies[i] = (double) i * div_frecuencies;

    powers[i] = 2.0 * div_power
                  * (data[2*i - 1]*data[2*i - 1] + data[2*i]*data[2*i]);
  }

  // For Nyquist frecuency
  i = Nx/2;
  frecuencies[i] = 1.0 / (2.0 * Delta);

  if( Nx % 2 == 0) {
    powers[i] = div_power * (data[Nx - 1] * data[Nx - 1]);
  }
  else {
    powers[i] = div_power
                  * (data[Nx - 2]*data[Nx - 2] + data[Nx - 1]*data[Nx - 1]);
  }
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

int main(int argc, char *argv[]) {
  /*============================================================
    Periodogram
  ============================================================

    The program uses the fast Fourier transform of GSL to
    calculate the periodogram of the time series.

    Output:

      periodogram.dat ---> file with the frecuencies vector
                            and its corresponding powers.

                          Output format:

                          column 1: frecuencies
                          columna 2: powers

    PROGRAM EXECUTION

      ./Periodogram.x Nx Delta

    Where:

      - Nx : number of sampling points.

      - Delta : Sampling interval time
  ============================================================*/

  // command line input variables
  int Nx;
  double Delta;

  // Program's internal variables
  int i;
  double *data=NULL, *frecuencies=NULL, *powers=NULL;

  // Output files
  FILE *file=NULL;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Nx = atoi(argv[1]);
  Delta  = atof(argv[2]);

  /*=========================================================
    Asign a memory block to the vectors:

    double *data;
    double *frecuencies;
    double *powers;

    Nx/2 is an integer division.
  =========================================================*/
  data  = (double *) malloc((size_t) Nx * sizeof(double));
  frecuencies = (double *) malloc((size_t) (Nx/2 + 1) * sizeof(double));
  powers  = (double *) malloc((size_t) (Nx/2 + 1) * sizeof(double));

  // Sampling points vector
  for(i = 0; i < Nx; i++) {
    data[i] = function(i*Delta);
  }

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
  gsl_fft_real_transform(data, 1, Nx, real, work);

  // Free memory
  gsl_fft_real_wavetable_free(real);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /*============================================================
    Periodogram
  ============================================================*/
  Periodogram(data,
              frecuencies,
              powers,
              Nx,
              Delta);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /*============================================================
    Print
  ============================================================*/
  file = fopen("periodogram.dat","w");

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
  free(data);
  free(frecuencies);
  free(powers);

return 0;
}
