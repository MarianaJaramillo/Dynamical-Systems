#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// Librería para manejo de matrices
#include "metodos_matrices.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  DEFINICIÓN DE FUNCIONES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

//============================================================

void System_Map(double *X,
                double *Y,
                double x,
                double y,
                double a,
                double b) {
  /*============================================================
    Retorna por referencia la solución al mapa de Hénon

    Hénon(x,y) = (a - x^2 - b y , x)
  ============================================================*/

  *X = a - x*x + b*y;
  *Y = x;
}

//============================================================

void Jacobian_Map(double **Jac,
                  double x,
                  double y,
                  double a,
                  double b) {
  /*============================================================
    Jacobian es la matriz (cuadrada) del jacobiano del
    mapa de Hénon:

      Jac[N][N] ---> donde N es el número de Variables

    Para el mapa de Hénon N=2

      Jac_Hénon(x,y) = [ [- 2 x, - b] , [1 , 0] ]

    En C:

      Jac[k][l] = *(*(Jac + k) + l)
  ============================================================*/
  //Jac[0][0] = - 2 * x
  *(*(Jac + 0) + 0) = - 2.0 * x;

  //Jac[0][1] = - b
  *(*(Jac + 0) + 1) = - b;

  //Jac[1][0] = 1
  *(*(Jac + 1) + 0) = 1.0;

  //Jac[1][1] = 0
  *(*(Jac + 1) + 1) = 0.0;
}


//============================================================

void Ortonormalizacion(double **basis,
                      double *normas_unnormalized_basis,
                      int dim) {
  /*============================================================
    Ortonormalizacion de la base almacenada en la matriz

      basis[dim][dim]

    Retorna por referencia la nueva base sobreescribiendo
    la matriz de la base original. También retorna las normas de
    los vectores de la base ortogonal sin normalizar
  ============================================================*/
  int i, j, k;
  double norma, producto;
  double *aux_vector1, *aux_vector2;

  /*============================================================
    Asignar bloques de memoria a los vectores

      *aux_vector1
      *aux_vector2
  ============================================================*/
  aux_vector1 = (double *) malloc((size_t) dim * sizeof(double));
  aux_vector2 = (double *) malloc((size_t) dim * sizeof(double));

  for(j = 0; j < dim; j++) {
    /*============================================================
      Copiar el vector j de la base en aux_vector1

        (Se copian todas las filas de la columna j)
    ============================================================*/
    for(i = 0; i < dim; i++) {
      aux_vector1[i] = basis[i][j];
    }

    /*============================================================
      Copiar los vectores

        k = j-1, j-2, ... 0

      de la base en aux_vector2

        (Se copian todas las filas de la columna k)

      Luego realiza el producto interno del vector j con cada
      vector k y realiza la resta vectorial
    ============================================================*/
    for(k = j-1; k > -1; k--) {
      // Copia el vector k
      for(i = 0; i < dim; i++) {
        aux_vector2[i] = basis[i][k];
      }

      // Producto interno
      producto = Producto_Interno(aux_vector1, aux_vector2, dim);

      // Resta vectorial del vector j con el vector k
      Restar_Vectores(1.0,
                      aux_vector1,
                      producto,
                      aux_vector2,
                      aux_vector1,
                      dim);
    }

    /*============================================================
      Calcular la norma del vector j y almacenar en el vector

        normas_unnormalized_basis
    ============================================================*/
    norma = Norma(aux_vector1, dim);

    normas_unnormalized_basis[j] = norma;

    /*============================================================
      Normalizo el vector j y lo ubico en la matriz de la base
    ============================================================*/
    for(i = 0; i < dim; i++) {
      basis[i][j] = aux_vector1[i] / norma;
    }
  }
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

int main(int argc, char *argv[]) {
  /*============================================================
    Lyapunov exponents
  ============================================================

    El programa calcula los exponentes de Lyapunov
    para el mapa de Henón.

    El programa genera el siguiente archivo de salida:

      Lyapunov_exponentes.dat ---> archivo con los exponentes
                                  de Lyapunov calculados hasta
                                  Nmax iteraciones.

    EJECUCION DEL PROGRAMA

      ./Lyapunov_exponents.x x0 y0 a b Nmax

    Donde:

      - x0 : condición inicial de x

      - y0 : condición inicial de y

      - a : valor de a

      - b : valor de b

      - Nmax : Número máximo de iteraciones
  ============================================================*/

  // Variables ingresadas por línea de comandos
  int Nmax;
  double x0, y0, a, b;

  // Variables internas del programa
  int i, iteracion, dim=2;
  double x, y, X, Y;
  double **Jac=NULL, **orthonormal_basis=NULL;
  double *normas_unnormalized_basis=NULL;
  double *Lyapunov_exponents=NULL;

  // Archivos de salida
  FILE *file=NULL;

  /*============================================================
    Abrir archivo para imprimir los datos
  ============================================================*/
  file = fopen("lyapunov.dat","w");

  fprintf(file,
    "#==================================================\n"
    "#  Datos generados para los exponentes de Lyapunov\n"
    "#==================================================\n"
    "#\n"
    "#   iteracion : columna 1 \n"
    "#           x : columna 2 \n"
    "#           y : columna 3 \n"
    "# exponente 1 : columna 4 \n"
    "# exponente 2 : columna 5 \n"
    "#\n"
  );

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  x0   = atof(argv[1]);
  y0   = atof(argv[2]);
  a    = atof(argv[3]);
  b    = atof(argv[4]);
  Nmax = atoi(argv[5]);

  /*============================================================
    Asignar bloque de memoria los vectores

      *normas_unnormalized_basis
      *Lyapunov_exponents
  ============================================================*/
  normas_unnormalized_basis = (double *) malloc((size_t) dim * sizeof(double));
  Lyapunov_exponents        = (double *) malloc((size_t) dim * sizeof(double));

  /*============================================================
    Asignar bloque de memoria a la matriz Jacobiana y a la base
    ortonormal

      **Jac
      **orthonormal_basis
  ============================================================*/
  Jac = (double **) malloc((size_t) dim * sizeof(double*));
  for(i = 0; i < dim; i++) {
    Jac[i] = (double *) malloc((size_t) dim * sizeof(double));
  }

  orthonormal_basis = (double **) malloc((size_t) dim * sizeof(double*));
  for(i = 0; i < dim; i++) {
    orthonormal_basis[i] = (double *) malloc((size_t) dim * sizeof(double));
  }

  /*============================================================
    Inicializo el vector Lyapunov_exponents
  ============================================================*/
  Lyapunov_exponents[0] = 0.0;
  Lyapunov_exponents[1] = 0.0;

  /*============================================================
    Inicializo la base ortonormal

    orthonormal_basis = ( 1 0
                          0 1 )
  ============================================================*/
  orthonormal_basis[0][0] = 1.0;
  orthonormal_basis[1][0] = 0.0;
  orthonormal_basis[0][1] = 0.0;
  orthonormal_basis[1][1] = 1.0;

  /*============================================================
    Inicializo 'x' y 'y'

    x = x0
    y = y0
  ============================================================*/
  x = x0;
  y = y0;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /*============================================================
    Ejecución
  ============================================================*/
  for(iteracion = 1; iteracion < (Nmax + 1); iteracion++) {
    // Calcular el jacobiano en el punto actual
    Jacobian_Map(Jac, x, y, a, b);

    // Aplicar el jacobiano a la base ortonormal
    Producto_Matrices(Jac, orthonormal_basis, dim, dim, dim, orthonormal_basis);

    // Ortonormalizar la nueva base
    Ortonormalizacion(orthonormal_basis, normas_unnormalized_basis, dim);

    for(i = 0; i < dim; i++) {
      Lyapunov_exponents[i] += log(normas_unnormalized_basis[i]);
    }

    /*============================================================
      Imprimir datos en disco
    ============================================================*/
    fprintf(file,"%d\t%lf\t%lf", iteracion, x, y);

    for(i = 0; i < dim; i++){
        fprintf(file,
          "\t%.16g"
          , Lyapunov_exponents[i] / (double) iteracion
        );
    }

    fprintf(file,"\n");

    // Evolucionar el sistema
    System_Map(&X, &Y, x, y, a, b);

    // Reiniciar condiciones actuales
    x = X;
    y = Y;

  }

  // Cerrar archivo
  fclose(file);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // Liberar memoria
  free(Jac);
  free(orthonormal_basis);
  free(normas_unnormalized_basis);
  free(Lyapunov_exponents);

return 0;
}
