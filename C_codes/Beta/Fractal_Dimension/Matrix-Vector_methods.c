#include "Matrix-Vector_methods.h"


//============================================================

void Vector_copy(double *a,
                double *a_asst,
                int n) {
  /*============================================================
    This function copies the elements from the original vector *a
    to the assistant vector *a_asst.
    (Optimized version using local memory)

    Both vectors must have the same size

      a[n]
      a_asst[n]

    In C:

      a[k] = *(a + k)
  ============================================================*/
  int i;

  for(i = 0; i < (n/3); i++) {
    *(a_asst + 3*i)     = *(a + 3*i);
    *(a_asst + 3*i + 1) = *(a + 3*i + 1);
    *(a_asst + 3*i + 2) = *(a + 3*i + 2);
  }

  for(i = 0; i < (n%3); i++) {
    *(a_asst + (n - 1) - i) = *(a + (n - 1) - i);
  }
}


//============================================================

void Matrix_copy(double **a,
                double **a_asst,
                int k,
                int l) {
  /*============================================================
    This function copies the elements from the original matrix **a
    to the assistant matrix **a_asst.
    (Optimized version using local memory)

    Both matrices must have the same dimensions

      a[k][l]
      a_asst[k][l]

    In C:

      a[k][l] = *(*(a + k) + l)
  ============================================================*/
  int i, j;

  for(i = 0; i < k; i++) {
    for(j = 0; j < (l/3); j++) {
      *(*(a_asst + i) + 3*j)     = *(*(a + i) + 3*j);
      *(*(a_asst + i) + 3*j + 1) = *(*(a + i) + 3*j + 1);
      *(*(a_asst + i) + 3*j + 2) = *(*(a + i) + 3*j + 2);
    }

    for(j = 0; j < (l%3); j++) {
      *(*(a_asst + i) + (l - 1) - j) = *(*(a + i) + (l - 1) - j);
    }
  }
}


//============================================================

void Vector_addition(double alpha,
                    double *a,
                    double beta,
                    double *b,
                    double *c,
                    int n) {
  /*============================================================
    This function executes the vectorial addition of the
    vectors 'a' and 'b', and returns the result by reference
    in the vector 'c'

      c = alpha * a + beta * b

    (Optimized version using local memory)

    The vectors must have the same size

      a[n]
      b[n]
      c[n]

    In C:

      a[k] = *(a + k)

    The work is done using assistant vectors for avoiding to
    overwrite the original vectors.
  ============================================================*/
  int i;
  double *a_asst=NULL, *b_asst=NULL;

  /*============================================================
    To assign memory blocks for the assistant vectors

      *a_asst
      *b_asst
  ============================================================*/
  a_asst = (double *) malloc((size_t) n * sizeof(double));
  b_asst = (double *) malloc((size_t) n * sizeof(double));

  /*============================================================
    To copy

      a[n] in a_asst[n]
      b[n] in b_asst[n]
  ============================================================*/
  Vector_copy(a, a_asst, n);
  Vector_copy(b, b_asst, n);

  for(i = 0; i < (n/3); i++) {
    *(c + 3*i)     =
                    *(a_asst + 3*i) * alpha
                    + *(b_asst + 3*i) * beta;
    *(c + 3*i + 1) =
                    *(a_asst + 3*i + 1) * alpha
                    + *(b_asst + 3*i + 1) * beta;
    *(c + 3*i + 2) =
                    *(a_asst + 3*i + 2) * alpha
                    + *(b_asst + 3*i + 2) * beta;
    }

  for(i = 0; i < (n%3); i++) {
    *(c + (n - 1) - i) =
                        *(a_asst + (n - 1) - i) * alpha
                        + *(b_asst + (n - 1) - i) * beta;
  }

  /*============================================================
    Free memory assigned to the assistant vectors

      *a_asst
      *b_asst
  ============================================================*/
  free(a_asst);
  free(b_asst);
}

//============================================================

void Vector_subtraction(double alpha,
                        double *a,
                        double beta,
                        double *b,
                        double *c,
                        int n) {
  /*============================================================
    This function executes the vectorial subtraction of the
    vectors 'a' and 'b' (a - b), and returns the result by reference
    in the vector 'c'

      c = alpha * a - beta * b

    (Optimized version using local memory)

    The vectors must have the same size

      a[n]
      b[n]
      c[n]

    In C:

      a[k] = *(a + k)

    The work is done using assistant vectors for avoiding to
    overwrite the original vectors.
  ============================================================*/
  int i;
  double *a_asst=NULL, *b_asst=NULL;

  /*============================================================
    To assign memory blocks for the assistant vectors

      *a_asst
      *b_asst
  ============================================================*/
  a_asst = (double *) malloc((size_t) n * sizeof(double));
  b_asst = (double *) malloc((size_t) n * sizeof(double));

  /*============================================================
    To copy

      a[n] in a_asst[n]
      b[n] in b_asst[n]
  ============================================================*/
  Vector_copy(a, a_asst, n);
  Vector_copy(b, b_asst, n);

  for(i = 0; i < (n/3); i++) {
    *(c + 3*i)     =
                    *(a_asst + 3*i) * alpha
                    - *(b_asst + 3*i) * beta;
    *(c + 3*i + 1) =
                    *(a_asst + 3*i + 1) * alpha
                    - *(b_asst + 3*i + 1) * beta;
    *(c + 3*i + 2) =
                    *(a_asst + 3*i + 2) * alpha
                    - *(b_asst + 3*i + 2) * beta;
  }

  for(i = 0; i < (n%3); i++) {
    *(c + (n - 1) - i) =
                        *(a_asst + (n - 1) - i) * alpha
                        - *(b_asst + (n - 1) - i) * beta;
  }

  /*============================================================
    Free memory assigned to the assistant vectors

      *a_asst
      *b_asst
  ============================================================*/
  free(a_asst);
  free(b_asst);
}

//============================================================

void Vector_Scalar_product(double *a,
                          double scalar,
                          int n,
                          double *c) {
  /*============================================================
    This function executes the product of the vector 'a' with
    a scalar and returns the result by reference in the vector 'c'

      c[n] = scalar * a[n]

    The work is done using assistant vectors for avoiding to
    overwrite the original ones.
  ============================================================*/
  double *a_asst=NULL;
  int k;

  /*============================================================
    To assign a memory block for the assistant vector

      *a_asst
  ============================================================*/
  a_asst = (double *) malloc((size_t) n * sizeof(double));

  /*============================================================
    To copy a[n] in a_asst[n]
  ============================================================*/
  Vector_copy(a, a_asst, n);

  /*============================================================
    To execute the product
  ============================================================*/
  for(k = 0; k < (n/3); k++) {
    c[3*k]     = scalar * a_asst[3*k];
    c[3*k + 1] = scalar * a_asst[3*k + 1];
    c[3*k + 2] = scalar * a_asst[3*k + 2];
  }

  for(k = 0; k < (n%3); k++) {
    c[(n - 1) - k] = scalar * a_asst[(n - 1) - k];
  }

  /*============================================================
    Free memory assigned to the assistant vector

      *a_asst
  ============================================================*/
  free(a_asst);
}


//============================================================

void Matrix_Scalar_product(double **a,
                          double scalar,
                          int m,
                          int n,
                          double **c) {
  /*============================================================
    This function executes the product of the matrix 'a' with
    a scalar and returns the result by reference in the
    matrix 'c'

      c[m][n] = scalar * a[m][n]


    The work is done using assistant matrices for avoiding to
    overwrite the original ones.
  ============================================================*/
  double **a_asst=NULL;
  int j, k;

  /*============================================================
    To assign a memory block for the assistant matrix

      **a_asst
  ============================================================*/
  a_asst = (double **) malloc((size_t) m * sizeof(double*));
  for(j = 0; j < m; j++) {
    a_asst[j] = (double *) malloc((size_t) n * sizeof(double));
  }

  /*============================================================
    To copy a[m][n] in a_asst[m][n]
  ============================================================*/
  Matrix_copy(a, a_asst, m, n);

  /*============================================================
    To execute the product
  ============================================================*/
  for(j = 0; j < m; j++) {
    // Optimized cycle in k
    for(k = 0; k < (n/3); k++) {
      c[j][3*k]     = scalar * a_asst[j][3*k];
      c[j][3*k + 1] = scalar * a_asst[j][3*k + 1];
      c[j][3*k + 2] = scalar * a_asst[j][3*k + 2];
    }

    for(k = 0; k < (n%3); k++) {
      c[j][(n - 1) - k] = scalar * a_asst[j][(n - 1) - k];
    }
  }

  /*============================================================
    Free memory assigned to the assistant matrix

      **a_asst
  ============================================================*/
  free(a_asst);
}


//============================================================

void Matrices_product(double **a,
                      double **b,
                      int m,
                      int l,
                      int n,
                      double **c) {
  /*============================================================
    This function executes the product of the matrices
    'a' and 'b' and returns the result by reference in the
    matrix 'c'

     a[m][l] * b[l][n] = c[m][n]

     The work is done using assistant matrices for avoiding to
     overwrite the original ones.
  ============================================================*/
  double add;
  double **a_asst=NULL, **b_asst=NULL;
  int i, j, k;

  /*============================================================
    To assign memory blocks for the assistant matrices

      **a_asst
      **b_asst
  ============================================================*/
  a_asst = (double **) malloc((size_t) m * sizeof(double*));
  for(i = 0; i < m; i++) {
    a_asst[i] = (double *) malloc((size_t) l * sizeof(double));
  }

  b_asst = (double **) malloc((size_t) l * sizeof(double*));
  for(i = 0; i < l; i++) {
    b_asst[i] = (double *) malloc((size_t) n * sizeof(double));
  }

  /*============================================================
    To copy a[m][l] in a_asst[m][l]
    To copy b[l][n] in b_asst[l][n]
  ============================================================*/
  Matrix_copy(a, a_asst, m, l);
  Matrix_copy(b, b_asst, l, n);

  /*============================================================
    To execute the product
  ============================================================*/
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++) {
      add = 0;
      // Optimized cycle in k
      for(k = 0; k < (l/3); k++) {
        add += a_asst[i][3*k]     * b_asst[3*k][j];
        add += a_asst[i][3*k + 1] * b_asst[3*k + 1][j];
        add += a_asst[i][3*k + 2] * b_asst[3*k + 2][j];
      }

      for(k = 0; k < (l%3); k++) {
        add += a_asst[i][(l - 1) - k] * b_asst[(l - 1) - k][j];
      }

      // To save the final result
      c[i][j] = add;
    }
  }

  /*============================================================
    Free memory assigned to the assistant matrices

      **a_asst
      **b_asst
  ============================================================*/
  free(a_asst);
  free(b_asst);
}


//============================================================

double Dot_product(double *a,
                  double *b,
                  int n) {
  /*============================================================
    This function executes the dot product (scalar product) of
    two vectors 'a' and 'b'

      a[n]
      b[n]

      dot product:   a . b = sum_i a[i]*b[i]
  ============================================================*/
  double product = 0.0;
  int i;

  /*============================================================
    To execute the product
  ============================================================*/
  for(i = 0; i < (n/3); i++) {
    product += a[3*i] * b[3*i];
    product += a[3*i + 1] * b[3*i + 1];
    product += a[3*i + 2] * b[3*i + 2];
  }

  for(i = 0; i < (n%3); i++) {
    product += a[(n - 1) - i] * b[(n - 1) - i];
  }

  return product;
}


//============================================================

double Norm(double *a,
            int n) {
  /*============================================================
    This function calculates the norm of the vector 'a'

      a[n]

      norm:   sqrt(a . a) = sqrt( sum_i a[i]*a[i] )
  ============================================================*/
  double product = 0.0;

  product = Dot_product(a, a, n);

  return sqrt(product);
}
