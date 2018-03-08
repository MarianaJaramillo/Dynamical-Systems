//
// Copyright (c) 2018 by Camilo-HG. All Rights Reserved.
//

// Library with defined structs
#include "Dynamical_Systems_Structs.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

int main(int argc, char *argv[]) {
  // Dynamical System Variables
  Dynamical_System sys;
  int dimension, Nparams, Npoints;
  //double initial_point[2], params[12];

  /*============================================================
    Setting out system's constant information
  ============================================================*/
  dimension = 2;
  Nparams = 12;
  transient = 1000000;
  Npoints = 1000000;

  Dynamical_System_alloc(&sys,
                        dimension,
                        Nparams,
                        transient,
                        Npoints);

  return 0;
}
