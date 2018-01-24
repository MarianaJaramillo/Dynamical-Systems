#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// Library for handling vectors and matrices
#include "Matrix-Vector_methods.h"

// Library with defined structs
#include "Dynamical_Systems_Structs.h"

/*************************************************************
  Methods and Routines
*************************************************************/

double logbase(double y, int b);

void Box_Dimension(Dynamical_System * sys,
                  double epsilon,
                  int max_exponent,
                  char *filename);
