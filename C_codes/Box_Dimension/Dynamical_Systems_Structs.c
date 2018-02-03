// Library with defined structs
#include "Dynamical_Systems_Structs.h"


//============================================================

void Dynamical_System_alloc(Dynamical_System * sys,
                            size_t dimension,
                            size_t Nparams,
                            size_t transient,
                            size_t Npoints) {
  // Method's internal variable
  int i;

  sys->dimension = dimension;
  sys->Nparams = Nparams;
  sys->Npoints = Npoints;
  sys->transient = transient;

  /*============================================================
    Setting up system's parameters
  ============================================================*/
  sys->params = (double *) malloc((size_t) sys->Nparams * sizeof(double));

  /*============================================================
    Setting up system's initial points
  ============================================================*/
  sys->initial_point = (double *) malloc((size_t) sys->dimension * sizeof(double));

  /*============================================================
    Setting up system's points in transient regime
  ============================================================*/
  sys->transient_points = (double **) malloc((size_t) sys->transient * sizeof(double*));
  for(i = 0; i < (sys->transient); i++) {
    sys->transient_points[i] = (double *) malloc((size_t) sys->dimension * sizeof(double));
  }

  /*============================================================
    Setting up system's points
  ============================================================*/
  sys->points = (double **) malloc((size_t) sys->Npoints * sizeof(double*));
  for(i = 0; i < (sys->Npoints); i++) {
    sys->points[i] = (double *) malloc((size_t) sys->dimension * sizeof(double));
  }
}


//============================================================

void Dynamical_System_free(Dynamical_System * sys) {
  /*============================================================
    Free system's initial points
  ============================================================*/
  free(sys->initial_point);

  /*============================================================
    Free system's parameters
  ============================================================*/
  free(sys->params);

  /*============================================================
    Free system's points in transient regime
  ============================================================*/
  free(sys->transient_points);

  /*============================================================
    Free system's points
  ============================================================*/
  free(sys->points);
}


//============================================================

void Dynamical_System_initialize(Dynamical_System * sys,
                                void * function,
                                const double * params,
                                const double * initial_point) {
  // Method's internal variable
  int i;

  /*============================================================
    Setting out system's function
  ============================================================*/
  sys->function = function;

  /*============================================================
    Setting out system's parameters
  ============================================================*/
  for(i = 0; i < sys->dimension; i++) {
    sys->params[i] = params[i];
  }

  /*============================================================
    Setting out system's initial points
  ============================================================*/
  for(i = 0; i < sys->dimension; i++) {
    sys->initial_point[i] = initial_point[i];
    sys->transient_points[0][i] = initial_point[i];
  }

  /*============================================================
    Setting out system's points in the transient regime
  ============================================================*/
  for(i = 1; i < (sys->transient); i++) {
    SYS_FN_EVOL(sys, sys->transient_points[i-1], sys->transient_points[i]);
  }

  /*============================================================
    Setting out system's points
  ============================================================*/
  for(i = 0; i < sys->dimension; i++) {
    sys->points[0][i] = sys->transient_points[sys->transient - 1][i];
  }

  for(i = 1; i < (sys->Npoints); i++) {
    SYS_FN_EVOL(sys, sys->points[i-1], sys->points[i]);
  }
}
