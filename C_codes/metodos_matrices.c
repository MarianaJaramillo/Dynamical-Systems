#include "metodos_matrices.h"


//============================================================

void Copiar_Vectores(double *a,
                    double *a_aux,
                    int n) {
  /*============================================================
    Copia los elementos de un vector a otro vector
    (versión optimizada aprovechando la memoria local)

    Ambos vectores deben ser del mismo tamaño

      a[n]
      a_aux[n]

    En C:

      a[k] = *(a + k)
  ============================================================*/
  int i;

  for(i = 0; i < (n/3); i++) {
    *(a_aux + 3*i)     = *(a + 3*i);
    *(a_aux + 3*i + 1) = *(a + 3*i + 1);
    *(a_aux + 3*i + 2) = *(a + 3*i + 2);
  }

  for(i = 0; i < (n%3); i++) {
    *(a_aux + (n - 1) - i) = *(a + (n - 1) - i);
  }
}


//============================================================

void Copiar_Matrices(double **a,
                    double **a_aux,
                    int k,
                    int l) {
  /*============================================================
    Copia los elementos de una matriz a otra matriz
    (versión optimizada aprovechando la memoria local)

    Ambas matrices deben ser del mismo tamaño

      a[k][l]
      a_aux[k][l]

    En C:

      a[k][l] = *(*(a + k) + l)
  ============================================================*/
  int i, j;

  for(i = 0; i < k; i++) {
    for(j = 0; j < (l/3); j++) {
      *(*(a_aux + i) + 3*j)     = *(*(a + i) + 3*j);
      *(*(a_aux + i) + 3*j + 1) = *(*(a + i) + 3*j + 1);
      *(*(a_aux + i) + 3*j + 2) = *(*(a + i) + 3*j + 2);
    }

    for(j = 0; j < (l%3); j++) {
      *(*(a_aux + i) + (l - 1) - j) = *(*(a + i) + (l - 1) - j);
    }
  }
}


//============================================================

void Sumar_Vectores(double alpha,
                    double *a,
                    double beta,
                    double *b,
                    double *c,
                    int n) {
  /*============================================================
    Realiza la suma vectorial de los vectores 'a' y 'b' y
    retorna el resultado por referencia en el vector 'c'

      c = alpha * a + beta * b

    Versión optimizada aprovechando la memoria local

    Los vectores deben ser del mismo tamaño

      a[n]
      b[n]
      c[n]

    En C:

      a[k] = *(a + k)

    Se trabaja con vectores auxiliares para evitar
    sobreescribir los vectores originales
  ============================================================*/
  int i;
  double *a_aux=NULL, *b_aux=NULL;

  /*============================================================
    Asignar bloque de memoria a los vectores auxiliar

      *a_aux
      *b_aux
  ============================================================*/
  a_aux = (double *) malloc((size_t) n * sizeof(double));
  b_aux = (double *) malloc((size_t) n * sizeof(double));

  /*============================================================
    Copiar

      a[n] en a_aux[n]
      b[n] en b_aux[n]
  ============================================================*/
  Copiar_Vectores(a, a_aux, n);
  Copiar_Vectores(b, b_aux, n);

  for(i = 0; i < (n/3); i++) {
    *(c + 3*i)     =
                    *(a_aux + 3*i) * alpha
                    + *(b_aux + 3*i) * beta;
    *(c + 3*i + 1) =
                    *(a_aux + 3*i + 1) * alpha
                    + *(b_aux + 3*i + 1) * beta;
    *(c + 3*i + 2) =
                    *(a_aux + 3*i + 2) * alpha
                    + *(b_aux + 3*i + 2) * beta;
    }

  for(i = 0; i < (n%3); i++) {
    *(c + (n - 1) - i) =
                        *(a_aux + (n - 1) - i) * alpha
                        + *(b_aux + (n - 1) - i) * beta;
  }

  /*============================================================
    Liberar memoria de los vectores auxiliares

      *a_aux
      *b_aux
  ============================================================*/
  free(a_aux);
  free(b_aux);
}

//============================================================

void Restar_Vectores(double alpha,
                    double *a,
                    double beta,
                    double *b,
                    double *c,
                    int n) {
  /*============================================================
    Realiza la resta vectorial de los vectores 'a' y 'b' (a - b)
    y retorna el resultado por referencia en el vector 'c'

      c = alpha * a - beta * b

    Versión optimizada aprovechando la memoria local

    Los vectores deben ser del mismo tamaño

      a[n]
      b[n]
      c[n]

    En C:

      a[k] = *(a + k)

    Se trabaja con vectores auxiliares para evitar
    sobreescribir los vectores originales
  ============================================================*/
  int i;
  double *a_aux=NULL, *b_aux=NULL;

  /*============================================================
    Asignar bloque de memoria a los vectores auxiliar

      *a_aux
      *b_aux
  ============================================================*/
  a_aux = (double *) malloc((size_t) n * sizeof(double));
  b_aux = (double *) malloc((size_t) n * sizeof(double));

  /*============================================================
    Copiar

      a[n] en a_aux[n]
      b[n] en b_aux[n]
  ============================================================*/
  Copiar_Vectores(a, a_aux, n);
  Copiar_Vectores(b, b_aux, n);

  for(i = 0; i < (n/3); i++) {
    *(c + 3*i)     =
                    *(a_aux + 3*i) * alpha
                    - *(b_aux + 3*i) * beta;
    *(c + 3*i + 1) =
                    *(a_aux + 3*i + 1) * alpha
                    - *(b_aux + 3*i + 1) * beta;
    *(c + 3*i + 2) =
                    *(a_aux + 3*i + 2) * alpha
                    - *(b_aux + 3*i + 2) * beta;
  }

  for(i = 0; i < (n%3); i++) {
    *(c + (n - 1) - i) =
                        *(a_aux + (n - 1) - i) * alpha
                        - *(b_aux + (n - 1) - i) * beta;
  }

  /*============================================================
    Liberar memoria de los vectores auxiliares

      *a_aux
      *b_aux
  ============================================================*/
  free(a_aux);
  free(b_aux);
}

//============================================================

void Producto_Vector_Escalar(double *a,
                            double escalar,
                            int n,
                            double *c) {
  /*============================================================
    Calcula el producto de un vector por un escalar

     c[n] = escalar * a[n]

     y retorna el resultado por referencia en el vector 'c'

    Se trabaja con vectores auxiliares para evitar
    sobreescribir los vectores originales
  ============================================================*/
  double *a_aux=NULL;
  int k;

  /*============================================================
    Asignar bloque de memoria a la matriz auxiliar

      *a_aux
  ============================================================*/
  a_aux = (double *) malloc((size_t) n * sizeof(double));

  /*============================================================
    Copiar a[n] en a_aux[n]
  ============================================================*/
  Copiar_Vectores(a, a_aux, n);

  /*============================================================
    Realizar el producto de vector por escalar
  ============================================================*/
  for(k = 0; k < (n/3); k++) {
    c[3*k]     = escalar * a_aux[3*k];
    c[3*k + 1] = escalar * a_aux[3*k + 1];
    c[3*k + 2] = escalar * a_aux[3*k + 2];
  }

  for(k = 0; k < (n%3); k++) {
    c[(n - 1) - k] = escalar * a_aux[(n - 1) - k];
  }

  /*============================================================
    Liberar memoria del vector auxiliar

      *a_aux
  ============================================================*/
  free(a_aux);
}


//============================================================

void Producto_Matriz_Escalar(double **a,
                            double escalar,
                            int m,
                            int n,
                            double **c) {
  /*============================================================
    Calcula el producto de una matriz por un escalar

     c[m][n] = escalar * a[m][n]

     y retorna el resultado por referencia en la matriz 'c'

    Se trabaja con matrices auxiliares para evitar
    sobreescribir las matrices originales
  ============================================================*/
  double **a_aux=NULL;
  int j, k;

  /*============================================================
    Asignar bloque de memoria a la matriz auxiliar

      **a_aux
  ============================================================*/
  a_aux = (double **) malloc((size_t) m * sizeof(double*));
  for(j = 0; j < m; j++) {
    a_aux[j] = (double *) malloc((size_t) n * sizeof(double));
  }

  /*============================================================
    Copiar a[m][n] en a_aux[m][n]
  ============================================================*/
  Copiar_Matrices(a, a_aux, m, n);

  /*============================================================
    Realizar el producto de matriz por escalar
  ============================================================*/
  for(j = 0; j < m; j++) {
    // Ciclo optimizado en k
    for(k = 0; k < (n/3); k++) {
      c[j][3*k]     = escalar * a_aux[j][3*k];
      c[j][3*k + 1] = escalar * a_aux[j][3*k + 1];
      c[j][3*k + 2] = escalar * a_aux[j][3*k + 2];
    }

    for(k = 0; k < (n%3); k++) {
      c[j][(n - 1) - k] = escalar * a_aux[j][(n - 1) - k];
    }
  }

  /*============================================================
    Liberar memoria de las matrices auxiliares

      **a_aux
  ============================================================*/
  free(a_aux);
}


//============================================================

void Producto_Matrices(double **a,
                      double **b,
                      int m,
                      int l,
                      int n,
                      double **c) {
  /*============================================================
    Calcula el producto de matrices

     a[m][l] * b[l][n] = c[m][n]

    Se trabaja con matrices auxiliares para evitar
    sobreescribir en las matrices originales
  ============================================================*/
  double sum;
  double **a_aux=NULL, **b_aux=NULL;
  int i, j, k;

  /*============================================================
    Asignar bloques de memoria a las matrices auxiliares

      **a_aux
      **b_aux
  ============================================================*/
  a_aux = (double **) malloc((size_t) m * sizeof(double*));
  for(i = 0; i < m; i++) {
    a_aux[i] = (double *) malloc((size_t) l * sizeof(double));
  }

  b_aux = (double **) malloc((size_t) l * sizeof(double*));
  for(i = 0; i < l; i++) {
    b_aux[i] = (double *) malloc((size_t) n * sizeof(double));
  }

  /*============================================================
    Copiar a[m][l] en a_aux[m][l]
    Copiar b[l][n] en b_aux[l][n]
  ============================================================*/
  Copiar_Matrices(a, a_aux, m, l);
  Copiar_Matrices(b, b_aux, l, n);

  /*============================================================
    Realizar el producto de matrices
  ============================================================*/
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++) {
      sum = 0;
      // Ciclo optimizado en k
      for(k = 0; k < (l/3); k++) {
        sum += a_aux[i][3*k]     * b_aux[3*k][j];
        sum += a_aux[i][3*k + 1] * b_aux[3*k + 1][j];
        sum += a_aux[i][3*k + 2] * b_aux[3*k + 2][j];
      }

      for(k = 0; k < (l%3); k++) {
        sum += a_aux[i][(l - 1) - k] * b_aux[(l - 1) - k][j];
      }

      // Almaceno el valor final
      c[i][j] = sum;
    }
  }

  /*============================================================
    Liberar memoria de las matrices auxiliares

      **a_aux
      **b_aux
  ============================================================*/
  free(a_aux);
  free(b_aux);
}


//============================================================

double Producto_Interno(double *a,
                      double *b,
                      int n) {
  /*============================================================
    Calcula el producto interno de los vectores a y b

      a[n]
      b[n]

      producto:   a . b = sum_i a[i]*b[i]
  ============================================================*/
  double producto = 0.0;
  int i;

  /*============================================================
    Realizar el producto
  ============================================================*/
  for(i = 0; i < (n/3); i++) {
    producto += a[3*i] * b[3*i];
    producto += a[3*i + 1] * b[3*i + 1];
    producto += a[3*i + 2] * b[3*i + 2];
  }

  for(i = 0; i < (n%3); i++) {
    producto += a[(n - 1) - i] * b[(n - 1) - i];
  }

  return producto;
}


//============================================================

double Norma(double *a,
            int n) {
  /*============================================================
    Calcula la norma del vector a

      a[n]

      norma:   sqrt(a . a) = sqrt( sum_i a[i]*a[i] )
  ============================================================*/
  double producto = 0.0;

  producto = Producto_Interno(a, a, n);

  return sqrt(producto);
}
