#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// Librerías GSL para hacer FFT
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  DEFINICIÓN DE FUNCIONES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

//============================================================

double function(double x) {
  double fx;

  fx = sin(x);

  return fx;
}

//============================================================

void Periodograma(double data[],
                  double frecuencias[],
                  double potencias[],
                  int Nx,
                  double Delta) {
  int i;

  double div_potencias, div_frecuencias;

  div_frecuencias = 1.0 / ( (double) Nx * Delta);
  div_potencias = 1.0/ (double) (Nx * Nx);

  // Para la frecuencia 0
  frecuencias[0] = 0.0;

  potencias[0] = div_potencias * (data[0]*data[0]);

  // Para las frecuencias internas
  for(i = 1; i < (Nx/2); i++){
    frecuencias[i] = (double) i * div_frecuencias;

    potencias[i] = 2.0 * div_potencias
                  * (data[2*i - 1]*data[2*i - 1] + data[2*i]*data[2*i]);
  }

  // Para la freuencia de Nyquist
  i = Nx/2;
  frecuencias[i] = 1.0 / (2.0 * Delta);

  if( Nx % 2 == 0) {
    potencias[i] = div_potencias * (data[Nx - 1] * data[Nx - 1]);
  }
  else {
    potencias[i] = div_potencias
                  * (data[Nx - 2]*data[Nx - 2] + data[Nx - 1]*data[Nx - 1]);
  }
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

int main(int argc, char *argv[]) {
  /*============================================================
    Periodograma
  ============================================================

    El programa utiliza la transformada discreta de Fourier
    de GSL (FTT) para calcular el periodograma de una serie
    de tiempo.

    El programa genera el siguiente archivo de salida:

      periodograma.dat ---> archivo con el vector de frecuencias
                            y sus correspondientes pesos
                            determinados por el periodograma.

    EJECUCION DEL PROGRAMA

      ./Periodograma.x Nx Delta

    Donde:

      - Nx : Número de puntos de puntos de muestreo

      - Delta : Intervalo de muestreo
  ============================================================*/

  // Variables ingresadas por línea de comandos
  int Nx;
  double Delta;

  // Variables internas del programa
  int i;
  double *data=NULL, *frecuencias=NULL, *potencias=NULL;

  // Archivos de salida
  FILE *file=NULL;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Nx = atoi(argv[1]);
  Delta  = atof(argv[2]);

  /*=========================================================
    Asignar un bloque de memoria a los arreglos:

    double *data;
    double *frecuencias;
    double *potencias;

    Nx/2 es una división entera
  =========================================================*/
  data  = (double *) malloc((size_t) Nx * sizeof(double));
  frecuencias = (double *) malloc((size_t) (Nx/2 + 1) * sizeof(double));
  potencias  = (double *) malloc((size_t) (Nx/2 + 1) * sizeof(double));

  // Arreglo con el muestreo de datos
  for(i = 0; i < Nx; i++) {
    data[i] = function(i*Delta);
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /*============================================================
    Variables para la FFT de GSL
  ============================================================*/
  gsl_fft_real_wavetable *real;
  gsl_fft_real_workspace *work;

  // Alocación de memoria
  work = gsl_fft_real_workspace_alloc(Nx);
  real = gsl_fft_real_wavetable_alloc(Nx);

  // FFT
  // Retorna data transformado
  gsl_fft_real_transform(data, 1, Nx, real, work);

  // Liberar memoria
  gsl_fft_real_wavetable_free(real);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /*============================================================
    Periodograma
  ============================================================*/
  Periodograma(data,
              frecuencias,
              potencias,
              Nx,
              Delta);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /*============================================================
    Imprimir datos en disco
  ============================================================*/
  file = fopen("periodograma.dat","w");

  fprintf(file,
    "#==================================================\n"
    "#  Datos generados por el periodograma\n"
    "#==================================================\n"
    "#\n"
    "# frecuencias : columna 1 \n"
    "#   potencias : columna 2 \n"
    "#\n"
  );

  for(i = 0; i < (Nx/2 + 1); i++){
      fprintf(file,
        "%.16g\t%.16g\n"
        , frecuencias[i], potencias[i]
      );
  }

  // Cerrar archivo
  fclose(file);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // Liberar memoria
  free(data);
  free(frecuencias);
  free(potencias);

return 0;
}
