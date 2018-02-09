//
// Copyright (c) 2018 by Camilo-HG. All Rights Reserved.
//

// Library for calculating fractal dimensions
#include "Fractal_Dimension.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FUNCTION DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

void System_Map(const double yo[],
                double f[],
                double *params) {
  /*============================================================
    Returns by reference the solution for the Hénon map:

      Hénon(x,y) = f(x,y) = (a - x^2 + b y , x)

    Hénon map doesn't depend on t.
  ============================================================*/
  double a, b;

  a = params[0];
  b = params[1];

  /*=========================================================
    Hénon map:
  =========================================================*/
  f[0] = a - yo[0]*yo[0] + b * yo[1];
  f[1] = yo[0];
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

int main(int argc, char *argv[]) {
  /*============================================================
    Henon map - Box dimension
  ============================================================

    This program calculates the fractal dimension for the Hénon
    map:

        Hénon(x,y) = f(x,y) = (a - x^2 + b y , x)

    OUTPUT FILES:

      Henon_box_dimension.dat ---> file which contains data
                                  generated by the program.

    PROGRAM EXECUTION

      ./Henon_Box_Dimension.x

    WHERE:

      - x0 : x's initial value

      - y0 : y's initial value

      - a : value of a

      - b : value of b

      - epsilon : grid's initial lenght

      - max_exponent : maximum reduction exponent
  ============================================================*/

  // Variables
  int max_exponent;
  int dimension, Nparams, transient, Npoints;
  double x0, y0, a, b, epsilon;
  double initial_point[2], params[2];
  Dynamical_System sys;
  double sys_Fractal_Dimension;
  double sys_Fractal_Dimension_SD;
  double tolerance_pct;
  int output_option, print_time_option;

  // Program's internal variables
  int i;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  x0 = 0.5;
  y0 = 0.5;
  a = 1.4;
  b = 0.3;

  /*============================================================
    Setting out system's constant information
  ============================================================*/
  dimension = 2;
  Nparams = 2;
  transient = 1000000;
  Npoints = 1000000;

  /*============================================================
    Alloc memory for Dynamical System
  ============================================================*/
  Dynamical_System_alloc(&sys,
                        (size_t) dimension,
                        (size_t) Nparams,
                        (size_t) transient,
                        (size_t) Npoints);

  /*============================================================
    Setting out system's initial point
  ============================================================*/
  initial_point[0] = x0;
  initial_point[1] = y0;

  /*============================================================
    Setting out system's initial point
  ============================================================*/
  params[0] = a;
  params[1] = b;

  /*============================================================
    Setting out configuration variables for Box_Dimension
  ============================================================*/
  epsilon  = 1.0;
  max_exponent = 13;
  tolerance_pct = 0.09;
  output_option = 1;
  print_time_option = 1;

  /*============================================================
    Initialize Dynamical System
  ============================================================*/
  Dynamical_System_initialize(&sys,
                              &System_Map,
                              params,
                              initial_point);

  /*============================================================
    Calculate Box_Dimension
  ============================================================*/
  Box_Dimension(&sys,
                &sys_Fractal_Dimension,
                &sys_Fractal_Dimension_SD,
                "Henon_Box_Dimension.dat",
                epsilon,
                max_exponent,
                tolerance_pct,
                output_option,
                print_time_option,
                "Henon_Box_Dimension_times_elapsed.dat");

  /*============================================================
    Save Dynamical System Dimension
  ============================================================*/
  FILE *file=NULL;
  file = fopen("Henon_Box_Dimension_Result.dat", "w");
  fprintf(file,
          "Henon map - Box Dimension = %.16g\n"
          "Henon map - Box Dimension standard deviation = %.16g\n"
          , sys_Fractal_Dimension, sys_Fractal_Dimension_SD
  );
  fclose(file);

  /*============================================================
    Save Dynamical System points
  ============================================================*/
  file = fopen("Henon_points.dat", "w");
  fprintf(file, "x\ty\n");
  for(i = 0; i < Npoints; i++) {
    fprintf(file,
            "%lf\t%lf\n"
            , sys.points[i][0], sys.points[i][1]
    );
  }
  fclose(file);

  /*============================================================
    Free memory for Dynamical System
  ============================================================*/
  Dynamical_System_free(&sys);

  return 0;
}
