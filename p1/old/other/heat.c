# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "mpi.h"

int main ( int argc, char *argv[] );
void heat_part ( int n, int p, int id, double x_min, double x_max );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
{
  double a = 0.0;
  double b = 1.0;
  int i;
  int id;
  int n;
  int p;
  double x_max;
  double x_min;
/*
  Startup:
*/ 
  MPI_Init ( &argc, &argv );

  MPI_Comm_rank ( MPI_COMM_WORLD, &id );

  MPI_Comm_size ( MPI_COMM_WORLD, &p );
/*
  Determine the portion of [A,B] to be assigned to processor ID.
*/ 
  n = 12;
  i = 0;

  x_min = ( ( double )( p * n + 1 - id * n - i ) * a   
          + ( double )(             id * n + i ) * b ) 
          / ( double ) ( p * n + 1              );

  i = n + 1;

  x_max = ( ( double )( p * n + 1 - id * n - i ) * a   
          + ( double )(             id * n + i ) * b ) 
          / ( double )( p * n + 1              );

  heat_part ( n, p, id, x_min, x_max );
/*
  Shut down.
*/ 
  MPI_Finalize ( );

  return p;
}
/******************************************************************************/

void heat_part ( int n, int p, int id, double x_min, double x_max )

/******************************************************************************/
/*
  HEAT_PART solves the heat equation over "part" of the domain.
*/ 
{
  double cfl;
  double *h;
  double *h_new;
  int i;
  int ierr;
  int j;
  int j_max;
  int j_min;
  double k;
  MPI_Status status;
  double t;
  double t_del;
  double t_max;
  double t_min;
  double t_new;
  int tag;
  double wtime;
  double *x;
  double x_del;

  h = ( double * ) malloc ( ( n + 2 ) * sizeof ( double ) );
  h_new = ( double * ) malloc ( ( n + 2 ) * sizeof ( double ) );
  x = ( double * ) malloc ( ( n + 2 ) * sizeof ( double ) );

  k = 0.002 / ( double ) p;
/*
  Time parameters:
*/ 
  j_min = 0;
  j_max = 100;
  t_min = 0.0;
  t_max = 10.0;
  t_del = ( t_max - t_min ) / ( double ) ( j_max - j_min );
/*
  Space parameters.
*/ 
  x_del = ( x_max - x_min ) / ( double ) ( n + 1 );
  for ( i = 0; i <= n + 1; i++ )
  {
    x[i] = ( ( double ) (         i ) * x_max   
           + ( double ) ( n + 1 - i ) * x_min ) 
           / ( double ) ( n + 1     );
  }
/*
  Set the initial value of H.
*/ 
  for ( i = 0; i <= n + 1; i++ )
  {
    h[i] = 95.0;
  }
/*
  Check the CFL condition.
*/ 
  cfl = k * t_del / x_del / x_del;

  if ( 0.5 <= cfl )
  {
    printf ( "  CFL condition failed.\n" );
    exit ( 1 );
  }
/*
  Each execution of this loop computes the solution at the next time.
*/ 
  wtime = MPI_Wtime ( );

  for ( j = 1; j <= j_max; j++ )
  {
/*
  Determine new time.
*/ 
    t_new = ( ( double ) (         j - j_min ) * t_max   
            + ( double ) ( j_max - j         ) * t_min ) 
            / ( double ) ( j_max     - j_min );
/*
  To set H_NEW(1:N), update the temperature based on the four point stencil.
*/ 
    for ( i = 1; i <= n; i++ )
    {
      h_new[i] = h[i] + t_del * ( 
        k * ( h[i-1] - 2.0 * h[i] + h[i+1] ) / x_del / x_del 
        + 2.0 * sin ( x[i] * t ) );
    }
/*
  MESSAGE #1:
    Processor ID sends its H_NEW(N) to neighbor ID+1...
*/ 
    tag = 1;

    if ( id < p - 1 )
    {
      MPI_Send ( &h_new[n], 1, MPI_DOUBLE, id+1, tag, MPI_COMM_WORLD );
    }
/*
  ...which receives it as H_NEW(0).
*/ 
    if ( 0 < id )
    {
      MPI_Recv ( &h_new[0], 1, MPI_DOUBLE, id-1, tag, MPI_COMM_WORLD, &status );
    }
    else
    {
      h_new[0] = 100.0 + 10.0 * sin ( t_new );
    }
/*
  Message #2:
    Processor ID sends its H_NEW(1) to neighbor ID-1...

  MISSING LINES OF CODE ARE INDICATED BY QUESTION MARKS!  INSERT THE CORRECT STATEMENTS.
*/ 
    tag = 2;

    if ( 0 < id )
    {
      ???
    }
/*
  ...which receives it as H_NEW(N+1).
*/ 
    if ( id < p - 1 )
    {
      ???
    }
    else
    {
      ???
    }
/*
  Update time and temperature.
*/ 
    t = t_new;

    for ( i = 0; i <= n + 1; i++ )
    {
      h[i] = h_new[i];
    }

  }

  wtime = MPI_Wtime ( ) - wtime;

  if ( id == 0 )
  {
    printf ( "\n" );
    printf ( "  Wall clock elapsed seconds = %f\n", wtime );
  }
/*
  Print solution for last value of T at end of computation.
*/  
  printf ( "%2d  T= %f\n", id, t );
  printf ( "%2d  X= ", id );
  for ( i = 0; i <= n + 1; i++ )
  {
    printf ( "%14f", x[i] );
  }
  printf ( "\n" );
  printf ( "%2d  H= ", id );
  for ( i = 0; i <= n + 1; i++ )
  {
    printf ( "%14f", h[i] );
  }
  printf ( "\n" );

  free ( h );
  free ( h_new );
  free ( x );

  return;
}
