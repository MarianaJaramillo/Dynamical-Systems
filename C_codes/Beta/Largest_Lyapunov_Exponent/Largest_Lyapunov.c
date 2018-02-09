// Library for calculating Largest Largest Lyapunov Exponent
#include "Largest_Lyapunov.h"


//============================================================

void Largest_Lyapunov_Exponent(const Dynamical_System * sys,
                              const double DRo,
                              double *Largest_Lyapunov) {
  int i;
  double Df_norm;
  double *R=NULL, *R1=NULL, *f=NULL, *f1=NULL, *Df=NULL;

  *Largest_Lyapunov = 0.0;

  R = (double *) malloc((size_t) sys->dimension * sizeof(double));
  R1 = (double *) malloc((size_t) sys->dimension * sizeof(double));
  f = (double *) malloc((size_t) sys->dimension * sizeof(double));
  f1 = (double *) malloc((size_t) sys->dimension * sizeof(double));
  Df = (double *) malloc((size_t) sys->dimension * sizeof(double));

  Vector_copy(sys->initial_point, R, sys->dimension);

  for(i = 0; i < sys->transient; i++) {
    SYS_FN_EVOL(sys, R, f);
    Vector_copy(f, R, sys->dimension);
  }

  Vector_Scalar_product(R,
                        1.0 + DRo,
                        sys->dimension,
                        R1);

  for(i = 0; i < (sys->Npoints); i++) {
    SYS_FN_EVOL(sys, R, f);
    SYS_FN_EVOL(sys, R1, f1);

    Vector_subtraction(1.0,
                      f1,
                      1.0,
                      f,
                      Df,
                      sys->dimension);

    Df_norm = Norm(Df,
                    sys->dimension);

    *Largest_Lyapunov += log(Df_norm/DRo);

    Vector_copy(f, R, sys->dimension);

    Vector_addition(1.0,
                    f,
                    DRo/Df_norm,
                    Df,
                    R1,
                    sys->dimension);
  }

  *Largest_Lyapunov = *Largest_Lyapunov/(double) (sys->Npoints);
}
