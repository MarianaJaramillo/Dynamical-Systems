// Library for calculating fractal dimensions
#include "Fractal_Dimension.h"


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

void Box_Dimension(Dynamical_System * sys,
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
