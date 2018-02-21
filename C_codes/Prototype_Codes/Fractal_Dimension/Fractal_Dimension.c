//
// Copyright (c) 2018 by Camilo-HG. All Rights Reserved.
//

// Library for calculating fractal dimensions
#include "Fractal_Dimension.h"


//============================================================

double logbase(double y, int b) {
  /*============================================================
    Returns the log base b of y
  ============================================================*/
  double lg;

  lg = log10(y) / log10(b);

  return lg;
}


//============================================================

void Box_Dimension(const Dynamical_System * sys,
                  double *sys_Fractal_Dimension,
                  double *sys_Fractal_Dimension_SD,
                  const char *filename,
                  const double epsilon,
                  const int max_exponent,
                  const double tolerance_pct,
                  const char *output_option,
                  const char *print_time_option,
                  const char *time_filename,
                  const char *method) {
  /*============================================================
    - *sys_Fractal_Dimension : pointer for System's Fractal Dimension.

    - *sys_Fractal_Dimension_SD : pointer for System's Fractal Dimension
                                  standard deviation.

    - filename: corresponds to output file name.

                Output format:

                  c1 = Decrease factor
                      Define the factor decrease according to the exponent
                      in the power of 2

                      factor_decrease = 2^exponent

                  c2 : Box linear side

                      box_linear_side = epsilon / factor_decrease

                  c3 : Nboxes

                  c4 : exponent

                  c5 : log2(1/box_linear_side)

                  c6 : log2(Nboxes)

                  c7 : contiguous slopes upwards

    - epsilon : grid's initial lenght

    - max_exponent : maximum reduction exponent
                    Exponents runs from 0 to max_exponent

    - output_option : 'false' = without Boxes_Data_Exponent_n.dat output (default)
                      'true' = with Boxes_Data_Exponent_n.dat output
                      (only works when sys->dimension == 2)

                      x0 =
                      x1 =
                      y0 =
                      y1 =

    - tolerance_pct : tolerance percentage for calculating
                      system's fractal dimension from contiguous
                      slopes (upwards)

    - print_time_option : 'false' = without Time_data.dat output (default)
                          'true' = with Time_data.dat

                          e = exponent
                          t = absolute time elapsed in seconds
                          h = hours elapsed
                          m = remained minutes elapsed
                          s = remained seconds elapsed

    - time_filename :

    - method : 'standard' (default)
               'CHmethod'
               'mixed'

  ============================================================*/
  double asst, l_asst, factor_decrease;
  int i, j, k, exponent, Npoints, Nboxes, newBoxes, point_detection;
  Hyperbox *Hyperboxes;
  double **box=NULL, **points=NULL;
  double *c1=NULL, *c2=NULL;
  int *c3=NULL, *c4=NULL;
  double *c5=NULL, *c6=NULL, *c7=NULL;
  double *selected_slopes=NULL;
  double diff_up, diff_down;
  int Nselected_slopes = 0;
  int any_coincidence = 0;
  double mean, sd;

  // time variables
  clock_t start;
  double execution_time, seconds, minutes, hours;
  FILE *time_file=NULL;

  /*============================================================
    Declare output file
  ============================================================*/
  FILE *file=NULL;
  file = fopen(filename,"w");

  /*============================================================
    Setting up output file

    c1 = Decrease factor
        Define the factor decrease according to the exponent
        in the power of 2

        factor_decrease = 2^exponent

    c2 : Box linear side

    c3 : Nboxes

    c4 : Exponent

    c5 : log2(1/box_linear_side)

    c6 : log2(Nboxes)

    c7 : contiguous slopes upwards
  ============================================================*/
  fprintf(file,
          "c1\tc2\tc3\tc4\tc5\tc6\tc7\n"
  );

  /*============================================================
    Declare boxfile
  ============================================================*/
  FILE *boxfile=NULL;
  char archivo[100];

  /*============================================================
    If print_time_option == 'true'

    Declare time variables

    Setting up time_file output

      e : exponent
      t : absolute time elapsed in seconds
      h : hours elapsed
      m : remained minutes elapsed
      s : remained seconds elapsed
  ============================================================*/
  switch(print_time_option) {
    case 'true' :
      // Open output file
      time_file = fopen(time_filename,"w");

      // Setting up time_file output
      fprintf(time_file,
              "e\tt\th\tm\ts\n");
      break;
    default :
      // print_time_option == 'false'
      break;
  }

  /*============================================================
    Asign a memory block for the assistant vectors

      *c1: for storing c1
      *c2: for storing c2
      *c3: for storing c3 - type int
      *c4: for storing c4 - type int
      *c5: for storing c5
      *c6: for storing c6
      *c7: for storing c7 ... Since it is only upwards, it has
                              an element less

      *selected_slopes : for storing selected slopes from c7
  ============================================================*/
  c1 = (double *) malloc((size_t) (max_exponent + 1) * sizeof(double));
  c2 = (double *) malloc((size_t) (max_exponent + 1) * sizeof(double));
  c3 = (int *) malloc((size_t) (max_exponent + 1) * sizeof(int));
  c4 = (int *) malloc((size_t) (max_exponent + 1) * sizeof(int));
  c5 = (double *) malloc((size_t) (max_exponent + 1) * sizeof(double));
  c6 = (double *) malloc((size_t) (max_exponent + 1) * sizeof(double));
  c7 = (double *) malloc((size_t) (max_exponent) * sizeof(double));
  selected_slopes = (double *) malloc((size_t) (max_exponent) * sizeof(double));

  /*============================================================
    Asign a memory block for the assistant matrix

      **points: for storing system's points and avoiding
                to modify it
  ============================================================*/
  points = (double **) malloc((size_t) (sys->Npoints) * sizeof(double*));
  for(j = 0; j < (sys->Npoints); j++) {
    points[j] = (double *) malloc((size_t) (sys->dimension) * sizeof(double));
  }

  /*============================================================
    Copy sys->points to **points
  ============================================================*/
  Matrix_copy(sys->points,
              points,
              sys->Npoints,
              sys->dimension);

  switch (method) {
    case 'mixed' :
      break;

    case 'CHmethod' :
      /*============================================================
        Asign a memory block for the assistant matrix

          **box: for storing the box's limits
      ============================================================*/
      box = (double **) malloc((size_t) (sys->dimension) * sizeof(double*));
      for(j = 0; j < (sys->dimension); j++) {
        box[j] = (double *) malloc((size_t) 2 * sizeof(double));
      }

      for(exponent = 0; exponent <= max_exponent; exponent++){
        /*============================================================
          Define the factor decrease according to the exponent in the
          power of 2

            factor_decrease = 2^exponent

          Remember:

            box_linear_side = epsilon / factor_decrease

          So:

            logbase_{2}(1 / box_linear_side) = logbase_{2}(2^exponent / epsilon)

          So, in units of epsilon:

            box_linear_side [epsilon] = 1.0 / factor_decrease

            logbase_{2}(1 / box_linear_side [epsilon])
              = logbase_{2}(1.0 / factor_decrease)
              = logbase_{2}(2^exponent) = exponentb

          And also:

            length / box_linear_side = factor_decrease * length / epsilon
        ============================================================*/
        factor_decrease = pow(2.0, (double) exponent);

        /*============================================================
          If print_time_option == 'true'
          Initialize execution time meassurement
        ============================================================*/
        switch(print_time_option) {
          case 'true' :
            start=clock();
            break;
          default :
            break;
        }

        /*============================================================
          Initialize Nboxes
        ============================================================*/
        Nboxes = 0;

        /*=========================================================
          If output_option == 'true'
          Setting up boxfile
        =========================================================*/
        switch(output_option) {
          case 'true':
            if(sys->dimension == 2) {
              sprintf(archivo,
                "Boxes_Data_Exponent_%d.dat"
                , exponent
              );
              boxfile = fopen(archivo,"w");

              fprintf(boxfile,
                      "x0\tx1\ty0\ty1\n");
            }
            break;
          default :
            // output_option == 'false'
            break;
        }

        /*============================================================
          Cycle for counting Nboxes
        ============================================================*/
        Npoints = sys->Npoints;
        while(Npoints > 0) {
          /*============================================================
            Increment number of boxes
          ============================================================*/
          Nboxes++;

          for(j = 0; j < sys->dimension; j++) {
            /*============================================================
              Take the first point and define box's limits
            ============================================================*/
            l_asst = factor_decrease * points[0][j] /epsilon;
            box[j][0] = floor(l_asst) * epsilon/factor_decrease;
            box[j][1] = epsilon/factor_decrease * (floor(l_asst) + 1.0);
          }

          /*=========================================================
            If output_option == 'true'
            Print in boxfile
          =========================================================*/
          switch(output_option) {
            case 'true':
              if(sys->dimension == 2) {
                fprintf(boxfile,
                        "%lf\t%lf\t%lf\t%lf\n"
                        , box[0][0], box[0][1], box[1][0], box[1][1]
                );
              }
              break;
            default :
              // output_option == 'false'
              break;
          }

          /*============================================================
            Obviously, the first point belongs to that box, so we move
            this point to the last position of the array of points
          ============================================================*/
          for(j = 0; j < sys->dimension; j++){
            asst = points[0][j];
            points[0][j] = points[Npoints - 1][j];
            points[Npoints - 1][j] = asst;
          }

          /*============================================================
            Decrease number of points
          ============================================================*/
          Npoints--;

          if (Npoints > 0) {
            for(i = 0; i < Npoints; i++) {
              /*============================================================
                Detect if the i-th point belongs to the box

                For every component j of the point, it must hold:

                  box[j][0] <= points[i][j] < box[j][1]

                The opposite is:

                  points[i][j] < box[j][0] OR box[j][1] <= points[i][j]
              ============================================================*/
              point_detection = 1;

              for(j = 0; j < sys->dimension; j++){
                if(points[i][j] < box[j][0] || box[j][1] <= points[i][j]){
                  point_detection = 0;
                  break;
                }
              }

              if(point_detection){
                /*============================================================
                  point_detection == 1

                  Then the point i belongs to the box. So we move
                  this point to the last position of the array of points
                ============================================================*/
                for(j = 0; j < sys->dimension; j++){
                  asst = points[i][j];
                  points[i][j] = points[Npoints - 1][j];
                  points[Npoints - 1][j] = asst;
                }

                /*============================================================
                  Decrease number of points
                ============================================================*/
                Npoints--;

                /*============================================================
                  Decrease counter i thus allowing to repeat evaluation for
                  the new i-th point
                ============================================================*/
                i--;
              }
            }
          }
        }
        /*=========================================================
          If output_option == 'true'
          Close boxfile
        =========================================================*/
        switch(output_option) {
          case 'true':
            if(sys->dimension == 2) {
              fclose(boxfile);
            }
            break;
          default :
            // output_option == 'false'
            break;
        }

        /*=========================================================
          Save data
        =========================================================*/
        c1[exponent] = factor_decrease;
        c2[exponent] = epsilon/factor_decrease;
        c3[exponent] = Nboxes;
        c4[exponent] = exponent;
        c5[exponent] = logbase(factor_decrease/epsilon, 2.0);
        c6[exponent] = logbase((double) Nboxes, 2.0);
        if(exponent > 0) {
          c7[exponent - 1] = (c6[exponent] - c6[exponent - 1])
                                          /(c5[exponent] - c5[exponent - 1]);
        }

        /*============================================================
          If print_time_option == 1
          Print in time_file the execution time

          e : exponent
          t : absolute time elapsed in seconds
          h : hours elapsed
          m : remained minutes elapsed
          s : remained seconds elapsed
        ============================================================*/
        switch(print_time_option) {
          case 'true' :
            execution_time = ((double) (clock() - start))/(double)CLOCKS_PER_SEC;
            hours = floor(execution_time/3600.);
            minutes = floor( ( (execution_time/3600.) - hours) * 60.0 );
            seconds = ( (execution_time/60.) - minutes ) * 60.0;

            fprintf(time_file,
                    "%d\t%.16g\t%lf\t%lf\t%lf\n"
                    , exponent, execution_time, hours, minutes, seconds
            );
          default :
            break;
        }
      }
      break;

    default :
      /*============================================================
        Initialize Nboxes
        ============================================================*/
        Nboxes = 0;

      /*============================================================
        Asign a memory block for the array of Hyperboxes

          *Hyperboxes: a vector of Hyperboxes for storing the boxes
                      and the points belonging to it
      ============================================================*/
      Hyperboxes = (Hyperbox *) malloc((size_t) 1 * sizeof(Hyperbox));

      /*============================================================
        factor_decrease = 2^exponent
      ============================================================*/
      factor_decrease = pow(2.0, (double) exponent);

      /*============================================================
        If print_time_option == 'true'
        Initialize execution time meassurement
      ============================================================*/
      switch(print_time_option) {
        case 'true' :
          start=clock();
          break;
        default :
          break;
      }

      /*=========================================================
        If output_option == 'true'
        Setting up boxfile
      =========================================================*/
      switch(output_option) {
        case 'true':
          if(sys->dimension == 2) {
            sprintf(archivo,
              "Boxes_Data_Exponent_%d.dat"
              , exponent
            );
            boxfile = fopen(archivo,"w");

            fprintf(boxfile,
                    "x0\tx1\ty0\ty1\n");
          }
          break;
        default :
          // output_option == 'false'
          break;
      }

      /*============================================================
        Cycle for counting Nboxes
        Made only for exponent == 0
      ============================================================*/
      exponent = 0;
      Npoints = sys->Npoints;
      while(Npoints > 0) {
        /*============================================================
          Increment number of boxes
        ============================================================*/
        Nboxes++;

        /*============================================================
          Set number of points in that box equal to 1
        ============================================================*/
        Hyperboxes[Nboxes - 1].Npoints = 1;

        /*============================================================
          Asign a memory block for the matrix Hyperboxes[].box

            *box: a vector of Hyperboxes for storing the boxes
                 and the points belonging to it
        ============================================================*/
        Hyperboxes[Nboxes - 1].box = (double **) malloc((size_t) (sys->dimension) * sizeof(double*));
        for(j = 0; j < (sys->dimension); j++) {
          Hyperboxes[Nboxes - 1].box[j] = (double *) malloc((size_t) 2 * sizeof(double));
        }

        /*============================================================
          Setting up box's points
        ============================================================*/
        Hyperboxes[Nboxes - 1].points =
          (double **) malloc((size_t) Hyperboxes[Nboxes - 1].Npoints * sizeof(double*));
        Hyperboxes[Nboxes - 1].points[Hyperboxes[Nboxes - 1].Npoints - 1] =
          (double *) malloc((size_t) sys->dimension * sizeof(double));

        /*============================================================
          Take the first point and define box's limits
        ============================================================*/
        for(j = 0; j < sys->dimension; j++) {
          l_asst = factor_decrease * points[0][j] /epsilon;
          Hyperboxes[Nboxes - 1].box[j][0] = floor(l_asst) * epsilon/factor_decrease;
          Hyperboxes[Nboxes - 1].box[j][1] = epsilon/factor_decrease * (floor(l_asst) + 1.0);
        }

        /*=========================================================
          If output_option == 'true'
          Print in boxfile
        =========================================================*/
        switch(output_option) {
          case 'true':
            if(sys->dimension == 2) {
              fprintf(boxfile,
                      "%lf\t%lf\t%lf\t%lf\n"
                      , Hyperboxes[Nboxes - 1].box[0][0], Hyperboxes[Nboxes - 1].box[0][1]
                      , Hyperboxes[Nboxes - 1].box[1][0], Hyperboxes[Nboxes - 1].box[1][1]
              );
            }
            break;
          default :
            // output_option == 'false'
            break;
        }

        /*============================================================
          Obviously, the first point belongs to that box, so we move
          this point to the last position of the array of points
          and we add this point in Hyperboxes[].point
        ============================================================*/
        for(j = 0; j < sys->dimension; j++){
          asst = points[0][j];
          points[0][j] = points[Npoints - 1][j];
          points[Npoints - 1][j] = asst;
          Hyperboxes[Nboxes - 1].points[Hyperboxes[Nboxes - 1].Npoints - 1][j] = asst;
        }

        /*============================================================
          Decrease number of points
        ============================================================*/
        Npoints--;

        if (Npoints > 0) {
          for(i = 0; i < Npoints; i++) {
            /*============================================================
              Detect if the i-th point belongs to the box

              For every component j of the point, it must hold:

                box[j][0] <= points[i][j] < box[j][1]

              The opposite is:

                points[i][j] < box[j][0] OR box[j][1] <= points[i][j]
            ============================================================*/
            point_detection = 1;

            for(j = 0; j < sys->dimension; j++){
              if(points[i][j] < boxes[Nboxes - 1][j][0]
                || boxes[Nboxes - 1][j][1] <= points[i][j]){
                point_detection = 0;
                break;
              }
            }

            if(point_detection){
              /*============================================================
                point_detection == 1

                Then the point i belongs to the box.

                Increment Hyperboxes.Npoints.

                Add this point to Hyperboxes[].points
              ============================================================*/
              Hyperboxes[Nboxes - 1].Npoints++;

              Hyperboxes[Nboxes - 1].points =
                (double **) realloc(Hyperboxes[Nboxes - 1].points, (size_t) Hyperboxes[Nboxes - 1].Npoints * sizeof(double*));
              Hyperboxes[Nboxes - 1].points[Hyperboxes[Nboxes - 1].Npoints - 1] =
                (double *) malloc((size_t) sys->dimension * sizeof(double));

              /*============================================================
                We move this point to the last position of the array
                of points
              ============================================================*/
              for(j = 0; j < sys->dimension; j++){
                asst = points[i][j];
                points[i][j] = points[Npoints - 1][j];
                points[Npoints - 1][j] = asst;
                Hyperboxes[Nboxes - 1].points[Hyperboxes[Nboxes - 1].Npoints - 1][j] = asst;
              }

              /*============================================================
                Decrease number of points
              ============================================================*/
              Npoints--;

              /*============================================================
                Decrease counter i thus allowing to repeat evaluation for
                the new i-th point
              ============================================================*/
              i--;
            }
          }
        }

        /*============================================================
          Realloc memory of array of boxes
        ============================================================*/
        boxes = (double ***) realloc((size_t) (Nboxes + 1) * sizeof(double**));
        boxes[Nboxes] = (double **) malloc((size_t) (sys->dimension) * sizeof(double*));
        for(j = 0; j < (sys->dimension); j++) {
          boxes[Nboxes][j] = (double *) malloc((size_t) 2 * sizeof(double));
        }
      }

      /*=========================================================
        If output_option == 'true'
        Close boxfile
      =========================================================*/
      switch(output_option) {
        case 'true':
          if(sys->dimension == 2) {
            fclose(boxfile);
          }
          break;
        default :
          // output_option == 'false'
          break;
      }

      /*=========================================================
        Save data
      =========================================================*/
      c1[exponent] = factor_decrease;
      c2[exponent] = epsilon/factor_decrease;
      c3[exponent] = Nboxes;
      c4[exponent] = exponent;
      c5[exponent] = logbase(factor_decrease/epsilon, 2.0);
      c6[exponent] = logbase((double) Nboxes, 2.0);
      if(exponent > 0) {
        c7[exponent - 1] = (c6[exponent] - c6[exponent - 1])
                                          /(c5[exponent] - c5[exponent - 1]);
      }

      /*============================================================
        If print_time_option == 1
        Print in time_file the execution time

        e : exponent
        t : absolute time elapsed in seconds
        h : hours elapsed
        m : remained minutes elapsed
        s : remained seconds elapsed
      ============================================================*/
      switch(print_time_option) {
        case 'true' :
          execution_time = ((double) (clock() - start))/(double)CLOCKS_PER_SEC;
          hours = floor(execution_time/3600.);
          minutes = floor( ( (execution_time/3600.) - hours) * 60.0 );
          seconds = ( (execution_time/60.) - minutes ) * 60.0;

          fprintf(time_file,
                  "%d\t%.16g\t%lf\t%lf\t%lf\n"
                  , exponent, execution_time, hours, minutes, seconds
          );
        default :
          break;
      }

      for(exponent = 1; exponent <= max_exponent; exponent++){
        for(k = 0; k < Nboxes; k++) {
          /*============================================================
            Realloc memory of array of boxes
            Remember we are dividing box linear side in halfs, so, we
            will have 2^(sys->dimension) new boxes, but as we already
            have the first box (the one that is divided), we only need
            to add 2^(sys->dimension) - 1 new boxes
          ============================================================*/
          newBoxes = (int) pow(2, sys->dimension) - 1;
          boxes = (double ***) realloc(
                            (size_t) (Nboxes + newBoxes) * sizeof(double**));
          for(i = 0; i < newBoxes; i++){
            boxes[Nboxes + i] = (double **) malloc((size_t) (sys->dimension) * sizeof(double*));
            for(j = 0; j < (sys->dimension); j++) {
              boxes[Nboxes + i][j] = (double *) malloc((size_t) 2 * sizeof(double));
            }
          }
        }
      }
      break;
  }

  /*=========================================================
    Print data in file
  =========================================================*/
  for(exponent = 0; exponent <= max_exponent; exponent++){
    if(exponent == 0) {
      fprintf(file,
            "%g\t%g\t%d\t%d\t%.16g\t%.16g\n"
            , c1[exponent], c2[exponent], c3[exponent], c4[exponent]
            , c5[exponent], c6[exponent]
      );
    }
    else {
      fprintf(file,
              "%g\t%g\t%d\t%d\t%.16g\t%.16g\t%.16g\n"
              , c1[exponent], c2[exponent], c3[exponent], c4[exponent]
              , c5[exponent], c6[exponent], c7[exponent - 1]
      );
    }
  }

  /*============================================================
    Calculate final_slope
    First two points and last point are eliminated
  ============================================================*/
  for(i = 2; i < (max_exponent - 1); i++) {

    diff_down = fabs(c7[i] - c7[i - 1])
              /c7[i];

    diff_up = fabs(c7[i] - c7[i + 1])
              /c7[i];

    if(diff_down < tolerance_pct && diff_up < tolerance_pct) {
      if(Nselected_slopes == 0) {
        any_coincidence = 1;
        Nselected_slopes += 2;
        selected_slopes[0] = c7[i - 1];
        selected_slopes[1] = c7[i];

      }
      else {
        selected_slopes[Nselected_slopes] = c7[i];
        Nselected_slopes++;
      }
    }
    else {
      if(any_coincidence) {
        any_coincidence = 0;
        selected_slopes[Nselected_slopes] = c7[i];
        Nselected_slopes++;
      }
    }
  }

  /*============================================================
    Calculate system's fractal dimension (and its standard
    deviation) as the mean of the selected slopes (and its
    standard deviation) using GSL library for statistics
  ============================================================*/
  mean = gsl_stats_mean(selected_slopes, 1, Nselected_slopes);
  sd = gsl_stats_sd_m(selected_slopes, 1, Nselected_slopes, mean);

  *sys_Fractal_Dimension = mean;
  *sys_Fractal_Dimension_SD = sd;

  /*============================================================
    Free **points
    Free **box
  ============================================================*/
  free(points);
  free(box);
  free(boxes);

  /*============================================================
    Free *c1
    Free *c2
    Free *c3
    Free *c4
    Free *c5
    Free *c6
    Free *c7
    Free *selected_slopes
  ============================================================*/
  free(c1);
  free(c2);
  free(c3);
  free(c4);
  free(c5);
  free(c6);
  free(c7);
  free(selected_slopes);

  /*============================================================
    Close file
  ============================================================*/
  fclose(file);

  /*============================================================
    If print_time_option == 1
    Close time_file
  ============================================================*/
  switch() {
    case 'true':
      fclose(time_file);
      break;
    default :
      // print_time_option == 'false'
      break;
  }
}
