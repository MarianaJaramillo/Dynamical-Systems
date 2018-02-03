#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

/*************************************************************
  Define structs
*************************************************************/

typedef struct Dynamical_System {
  /*============================================================
    STRUCT: Dynamical_System

      System[i] = function(System[i-1])

        function
        ========
        void (*function) (const double point[],
                          double f[],
                          const double *params);

        params[Nparams]
        ==========================


        points[Npoints][dimension]
        ==========================

          **points: Stores system's time series.

          Npoints: Gives us information about how
                  many points has the time series.

          dimension: Gives us information about the
                    dimensionality of the system.

          Thus:

            points[i][j] ---> corresponds to the j-th component
                              of the i-th point.


        initial_point[dimension]
        ==========================

          *initial_point: Stores initial point of the system.


        transient_points[transient][dimension]
        ======================================

          **transient_points: Stores the transient points of the system.

          transient: Gives us information about how
                    many points has the time series in the transient regime.

          Thus:

            transient_points[i][j] ---> corresponds to the j-th component
                                        of the i-th point in the transient
                                        regime.

          Note that:

            transient_points[0][j] = initial_point[j]
            points[0][j] = transient_points[transient - 1][j]
  ============================================================*/
  void (*function) (const double initial_point[],
                    double f[],
                    const double *params);
  double *params;
  double *initial_point;
  double **transient_points;
  double **points;
  size_t dimension;
  size_t Nparams;
  size_t transient;
  size_t Npoints;
} Dynamical_System;

/*************************************************************
  Methods and Routines
*************************************************************/

/* Function evaluation macros */
#define SYS_FN_EVOL(S, yo, f)  (*((S)->function))(yo,f,(S)->params)

void Dynamical_System_alloc(Dynamical_System * sys,
                            size_t dimension,
                            size_t Nparams,
                            size_t transient,
                            size_t Npoints);

void Dynamical_System_free(Dynamical_System * sys);

void Dynamical_System_initialize(Dynamical_System * sys,
                                void * function,
                                const double * params,
                                const double * initial_point);
