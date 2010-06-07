#include "mpi.h" 
#include <stdio.h> 
#include <math.h> 

#define TRUE 1
int
main( int argc, char *argv[] ) 
{ 
    int n, myid, numprocs, i; 
    double PI25DT = 3.141592653589793238462643; 
    double mypi, pi, h, sum, x, start, end, time, glob_time; 

    MPI_Init(&argc,&argv); 
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
    MPI_Comm_rank(MPI_COMM_WORLD,&myid); 

    while (1) { 
        if (myid == 0) { 
            printf("enter the number of intervals: (0 quits) "); 
            scanf("%d",&n); 
        } 
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); 
        if (n == 0) 
            break; 
        else { 
            h   = 1.0 / (double) n; 
            sum = 0.0;
	    start = MPI_Wtime();
            for (i = myid + 1; i <= n; i += numprocs) { 
                x = h * ((double)i - 0.5); 
                sum += (4.0 / (1.0 + x*x)); 
            }
            mypi = h * sum;
            end = MPI_Wtime();
            time += end - start;
	    MPI_Reduce(&time, &glob_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	    MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
            if (myid == 0)  
                printf("pi is approx %.16f, error is %.16f and time is %.16f\n", pi, fabs(pi - PI25DT),end-start); 
        } 
    } 
    MPI_Finalize(); 
    return 0; 
} 
