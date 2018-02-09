#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// Library for handling vectors and matrices
#include "Matrix-Vector_methods.h"

// Library with defined structs
#include "Dynamical_Systems_Structs.h"

#ifndef LYAPUNOV_EXPONENT  // header guard. Unique identifier for each header
#define LYAPUNOV_EXPONENT

/*************************************************************
  Methods and Routines
*************************************************************/

void Largest_Lyapunov_Exponent(const Dynamical_System * sys,
                              const double DRo,
                              double *Largest_Lyapunov);

#endif // LYAPUNOV_EXPONENT
