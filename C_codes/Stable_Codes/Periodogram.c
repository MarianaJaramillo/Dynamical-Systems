//
// Copyright (c) 2018 by Camilo-HG. All Rights Reserved.
//

// Library for handling vectors and matrices
#include "Periodogram.h"


//============================================================
void DynamicalSystem2array(Dynamical_System * sys,
                          const int coordinate,
                          double *data,
                          int *Ndata) {
  /*============================================================
    Transforms data from Dynamical_System' points to data array

    - *data = NULL
              data's vector

    - Ndata : number of sampling points.

    Remember: sys.points[Npoints][dimension]
  ============================================================*/
  int i;

  *Ndata = sys->Npoints;

  data = (double *) malloc((size_t) (*Ndata) * sizeof(double));

  for(i = 0; i < (*Ndata); i++) {
    data[i] = sys->points[i][coordinate];
  }
}


//============================================================

void Periodogram(const double data[],
                const int Ndata,
                const double Delta,
                double *frecuencies,
                double *powers,
                int *Nfrecuencies) {
  /*============================================================
    The program uses the fast Fourier transform of GSL to
    calculate the periodogram of the time series.

    - data[]: data's vector

    - Ndata : number of sampling points.

    - Delta : Sampling interval time

    - *frecuencies = NULL

    - *powers = NULL

    - *Nfrecuencies = number of frecuencies returned
  ============================================================*/

  // Program's internal variables
  int i;
  double *data_aux=NULL;
  double div_power, div_frecuencies;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /*=========================================================
    Nfrecuencies = Ndata/2 + 1

    Asign a memory block to the vectors:

      double *data_aux;
      double *frecuencies;
      double *powers;
  =========================================================*/
  *Nfrecuencies = Ndata/2 + 1;

  data_aux = (double *) malloc((size_t) Ndata * sizeof(double));
  frecuencies = (double *) malloc((size_t) Nfrecuencies * sizeof(double));
  powers  = (double *) malloc((size_t) Nfrecuencies * sizeof(double));

  /*=========================================================
    Copy data to data_aux to avoid to modify data
  =========================================================*/
  Vector_copy(data,
              data_aux,
              Ndata);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /*============================================================
    FFT variables
  ============================================================*/
  gsl_fft_real_wavetable *real;
  gsl_fft_real_workspace *work;

  // memory allocation
  work = gsl_fft_real_workspace_alloc(Ndata);
  real = gsl_fft_real_wavetable_alloc(Ndata);

  // FFT
  // Returns transformed data
  gsl_fft_real_transform(data_aux, 1, Ndata, real, work);

  // Free memory
  gsl_fft_real_wavetable_free(real);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /*============================================================
    Periodogram
  ============================================================*/

  div_frecuencies = 1.0 / ( (double) Ndata * Delta);
  div_power = 1.0/ (double) (Ndata * Ndata);

  // For frecuency 0
  frecuencies[0] = 0.0;

  powers[0] = div_power * (data_aux[0]*data_aux[0]);

  // For internal frecuencies
  for(i = 1; i < (Ndata/2); i++){
    frecuencies[i] = (double) i * div_frecuencies;

    powers[i] = 2.0 * div_power
                  * (data_aux[2*i - 1]*data_aux[2*i - 1] + data_aux[2*i]*data_aux[2*i]);
  }

  // For Nyquist frecuency
  i = Ndata/2;
  frecuencies[i] = 1.0 / (2.0 * Delta);

  if( Ndata % 2 == 0) {
    powers[i] = div_power * (data_aux[Ndata - 1] * data_aux[Ndata - 1]);
  }
  else {
    powers[i] = div_power
                  * (data_aux[Ndata - 2]*data_aux[Ndata - 2] + data_aux[Ndata - 1]*data_aux[Ndata - 1]);
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // Free memory
  free(data_aux);
}
