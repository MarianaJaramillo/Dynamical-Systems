#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// Librería para manejo de matrices
#include "metodos_matrices.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  DEFINICIÓN DE FUNCIONES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*============================================================
  System's Time Series

    points[Npoints][dimension]

      **points: It stores the points of the system.

      Npoints: It gives us information about how
              many points has the time series.

      dimension: It gives us information about the
                dimensionality of the system.

    Thus:

      points[i][j] ---> corresponds to the j-th component
                        of the i-th point.
*/
typedef struct Time_Series {
  double **points;
  size_t Npoints;
  size_t dimension;
} Time_Series;


//============================================================

double logbase(double y, int b) {
  /*============================================================
    Returns the log base b of y
  ============================================================*/
  double lg;

  lg = log10(y) / log10(b);

  return lg;
}


//============================================================

void Box_Dimension(Time_Series * sys,
                  double epsilon,
                  int max_exponent,
                  char *filename) {
  /*============================================================
    - filename: corresponds to output file name.

    - max_exponent:
    - epsilon:
  ============================================================*/
  double aux, l_aux, factor_decrease;
  int i, j, exponent, Npuntos, Nboxes, point_detection;
  double **box=NULL;
  FILE *file=NULL;

  file = fopen(filename,"w");

  fprintf(file,
    "#==================================================\n"
    "#  Box Dimension - Generated Data \n"
    "#==================================================\n"
    "#\n"
    "# Epsilon used : %g \n"
    "# Maximum exponent = %d \n"
    "#\n"
    "# Remember: \n"
    "#\n"
    "#    box_linear_side = epsilon / factor_decrease \n"
    "#\n"
    "# So:\n"
    "#\n"
    "#    log2(1 / box_linear_side) = log2(2^exponent / epsilon) \n"
    "#\n"
    "# So, in units of epsilon: \n"
    "# \n"
    "#    box_linear_side [epsilon] = 1.0 / factor_decrease \n"
    "#\n"
    "#    logbase_{2}(1 / box_linear_side [epsilon]) \n"
    "#      = logbase_{2}(1.0 / factor_decrease) \n"
    "#      = logbase_{2}(2^exponent) = exponent \n"
    "#\n"
    "# And also:\n"
    "#\n"
    "#  length / box_linear_side = factor_decrease * length / epsilon\n"
    "#\n"
    "#         Decrease factor : column 1 \n"
    "#         Box linear side : column 2 \n"
    "#                  Nboxes : column 3 \n"
    "#                Exponent : column 4 \n"
    "# log2(1/box_linear_side) : column 5 \n"
    "#            log2(Nboxes) : column 6 \n"
    "#\n"
    , epsilon, max_exponent
  );

  /*============================================================
    Asign a memory block for the auxiliar matrix

      **box

    which will store the box's limits
  ============================================================*/
  box = (double **) malloc((size_t) (sys->dimension) * sizeof(double*));
  for(j = 0; j < (sys->dimension); j++) {
    box[j] = (double *) malloc((size_t) 2 * sizeof(double));
  }

  for(exponent = 0; exponent <= max_exponent; exponent++){
    /*============================================================
      Define the factor decrease according to the exponent in the
      power of 2

        factor_decrease = 2^exponent

      Remember:

        box_linear_side = epsilon / factor_decrease

      So:

        logbase_{2}(1 / box_linear_side) = logbase_{2}(2^exponent / epsilon)

      So, in units of epsilon:

        box_linear_side [epsilon] = 1.0 / factor_decrease

        logbase_{2}(1 / box_linear_side [epsilon])
          = logbase_{2}(1.0 / factor_decrease)
          = logbase_{2}(2^exponent) = exponent

      And also:

        length / box_linear_side = factor_decrease * length / epsilon
    ============================================================*/
    factor_decrease = pow(2.0, (double) exponent);

    /*============================================================
      Initialize Nboxes
    ============================================================*/
    Nboxes = 0;

    /*============================================================
      Cycle for counting Nboxes
    ============================================================*/
    Npuntos = sys->Npoints;
    while(Npuntos > 0) {
      for(j = 0; j < sys->dimension; j++) {
        /*============================================================
          Take the first point and define box's limits
        ============================================================*/
        l_aux = factor_decrease * sys->points[0][j] /epsilon;
        box[j][0] = floor(l_aux) * epsilon/factor_decrease;
        box[j][1] = epsilon/factor_decrease * (floor(l_aux) + 1.0);
      }
      /*============================================================
        Increment number of boxes
      ============================================================*/
      Nboxes++;

      /*============================================================
        Obviously, the first point belongs to that box, so we move
        this point to the last position of the array of points
      ============================================================*/
      for(j = 0; j < sys->dimension; j++){
        aux = sys->points[0][j];
        sys->points[0][j] = sys->points[Npuntos - 1][j];
        sys->points[Npuntos - 1][j] = aux;
      }

      /*============================================================
        Decrease number of points
      ============================================================*/
      Npuntos--;

      if (Npuntos > 0) {
        for(i = 0; i < Npuntos; i++) {
          /*============================================================
            Detect if the i-th point belongs to the box

            For every component j of the point, it must hold:

              box[j][0] <= points[i][j] < box[j][1]

            The opposite is:

              points[i][j] < box[j][0] OR box[j][1] <= points[i][j]
          ============================================================*/
          point_detection = 1;

          for(j = 0; j < sys->dimension; j++){
            if(sys->points[i][j] < box[j][0] || box[j][1] <= sys->points[i][j]){
              point_detection = 0;
              break;
            }
          }

          if(point_detection){
            /*============================================================
              point_detection == 1

              Then the point i belongs to the box. So we move
              this point to the last position of the array of points
            ============================================================*/
            for(j = 0; j < sys->dimension; j++){
              aux = sys->points[i][j];
              sys->points[i][j] = sys->points[Npuntos - 1][j];
              sys->points[Npuntos - 1][j] = aux;
            }

            /*============================================================
              Decrease number of points
            ============================================================*/
            Npuntos--;

            /*============================================================
              Decrease counter i thus allowing to repeat evaluation for
              the new i-th point
            ============================================================*/
            i--;
          }
        }
      }
    }

    fprintf(file,
            "%g\t%g\t%d\t%d\t%.16g\t%.16g\n"
            , factor_decrease, epsilon/factor_decrease, Nboxes, exponent
            , logbase(factor_decrease/epsilon, 2.0)
            , logbase((double) Nboxes, 2.0));
  }
  fclose(file);
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

int main(int argc, char *argv[]) {
  /*============================================================
    Largest Lyapunov exponent
  ============================================================

    El programa calcula la dimensión fractal por conteo de cajas
    para el mapa de Henón.

    El programa genera el siguiente archivo de salida:

      Henon_box_dimension.dat ---> archivo con el exponente más largo
                                  de Lyapunov calculado hasta
                                  Nmax iteraciones.

    EJECUCION DEL PROGRAMA

      ./Box_Dimension.x x0 y0 a b Nmax DRo

    Donde:

      - x0 : condición inicial de x

      - y0 : condición inicial de y

      - a : valor de a

      - b : valor de b

      - Nmax : Número máximo de iteraciones

      - epsilon : tamaño base de la grilla

      - max_exponent : exponente máximo de la reducción
  ============================================================*/

  // Variables
  int max_exponent;
  double x0, y0, a, b, epsilon;
  Time_Series sys;

  // Variables internas del programa
  int i;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  x0 = 0.5;
  y0 = 0.5;
  a = 1.4;
  b = 0.3;
  epsilon  = 0.1;
  max_exponent = 10;

  sys.dimension = 2;
  sys.Npoints = 10000;
  sys.points = (double **) malloc((size_t) sys.Npoints * sizeof(double*));
  for(i = 0; i < (sys.Npoints); i++) {
    sys.points[i] = (double *) malloc((size_t) sys.dimension * sizeof(double));
  }

  sys.points[0][0] = x0;
  sys.points[0][1] = y0;

  for(i = 1; i < (sys.Npoints); i++) {
    sys.points[i][0] = a - sys.points[i-1][0]*sys.points[i-1][0]
                      + b * sys.points[i-1][1];
    sys.points[i][1] = sys.points[i-1][0];
  }

  Box_Dimension(&sys,
                epsilon,
                max_exponent,
                "Henon_box_dimension.dat");

  return 0;
}
