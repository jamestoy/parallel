/* Sequential Solution to Steady-State Heat Problem */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>

#define TAG 0 

int
main (int argc, char *argv[])
{
    int n; double eps;
    if(argc != 3) {
        printf("usage: ssheat <grid size> <epsilon>\n");
        return -1;
    } else {
        n = atoi(argv[1]);
        eps = atof(argv[2]);
    }

    double diff,
           tdiff,
           global_diff;    /* Change in value */
    
    int i,j,
        iter,
        lr,hr,local,
        id,p,
        start,end;

    double mean;                 /* Average boundary value */
    
    double **u;
    double **w;
    double *_n,
           *_s;

    MPI_Status stat;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    /* Method described in class discussion
     * uses:
     * 1) allocation step
     * 2) distinguishing start / final in loops
     */
    lr = id * n / p;
    hr = (id + 1) * n / p - 1;
    local = hr - lr + 1;

    /* Allocation step */
    /*assert((u = malloc(local * sizeof(double *))) != NULL);
    assert((w = malloc(local * sizeof(double *))) != NULL);
    for(i = 0; i < local; i++) { //XXX:(jft) s/n/local
        assert((u[i] = malloc(n * sizeof(double))) != NULL);
        assert((w[i] = malloc(n * sizeof(double))) != NULL);
    }*/
    /* Ghost Rows */
    assert((_n = malloc(n * sizeof(double))) != NULL);
    assert((_s = malloc(n * sizeof(double))) != NULL);

    double *ualloc,
           *walloc;
    assert((ualloc = (double *) malloc(local * n * sizeof(double))) != NULL);
    assert((walloc = (double *) malloc(local * n * sizeof(double))) != NULL);
    assert((u = (double **) malloc(local * sizeof(double *))) != NULL);
    assert((w = (double **) malloc(local * sizeof(double *))) != NULL);
    for (i = 0; i < local; i++) {
        w[i] = &walloc[i*n];
        u[i] = &ualloc[i*n];
    }

    double t = MPI_Wtime();
    
    /* Set boundary values and compute mean boundary value */
    mean = 0.0;
    for (i = 0; i < n; i++) {
        u[i][0] = u[i][n-1] = u[0][i] = 100.0;
        u[n-1][i] = 0.0;
        mean += u[i][0] + u[i][n-1] + u[0][i] + u[n-1][i];
    }
    mean /= (4.0 * n);

    /* Initialize interior values */
    for(i = 1; i < n-1; i++)
        for(j = 1; j < n-1; j++) u[i][j] = mean;

    /* Set iterator boundaries */
    if(lr == 0)
        start = 1;          
    else 
        start = 0;          
    if(hr == n-1) 
        end = local-2;
    else
        end = local-1;

    /* Compute the steady-state solution */
    iter = 0;
    for(;;) {
        diff = 0.0;
        
        if (id != p--)
            MPI_Send(w[local--], n, MPI_DOUBLE, id+1, TAG, MPI_COMM_WORLD);
        if (id != 0) {
            MPI_Recv(_n, n, MPI_DOUBLE, id-1, TAG, MPI_COMM_WORLD, &stat);
            MPI_Send(w[0], n, MPI_DOUBLE, id-1, TAG, MPI_COMM_WORLD);
        }
        if (id != p-1)
            MPI_Recv(_s, n, MPI_DOUBLE, id+1, TAG, MPI_COMM_WORLD, &stat);
#pragma omp parallel private(i,j,tdiff)
{
        tdiff = 0.0;
#pragma omp for
        for(i = start; i <= end; i++)    //XXX:(jft) s/</<=/g
            for(j = 1; j < n-1; j++) {
                if(i == 0)
                    w[i][j] = (_n[j] + u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]) / 4.0;
                else if(i == local - 1)
                    w[i][j] = (_s[j] + u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]) / 4.0;
                else                     // regular computation
                    w[i][j] = (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]) / 4.0;
                
                /* Delimeter */
                if (fabs(w[i][j] - u[i][j]) > tdiff) //XXX: (jft) s/diff/tdiff
                    tdiff = fabs(w[i][j] - u[i][j]);
            }
#pragma omp for nowait
        for(i = 1; i < n-1; i++)
            for(j = 1; j < n-1; j++)
                u[i][j] = w[i][j];
#pragma omp critical
        if (tdiff > diff) diff = tdiff;
}
        MPI_Allreduce(&diff, &global_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (global_diff <= eps) break;
        iter++;
    }
    /* Print number of iteration */
    printf("number of iterations: %d\n\n",iter);
    
    /* Derive run-time */
    t += MPI_Wtime();

    /* Print solution */
    /*for (i = 0; i < n; i++) {
        for(j = 0; j < n; j++)
            printf("%6.2f ", u[i][j]);
        putchar('\n');
    }*/
    
    //free(u); free(w); free(_n); free(_s);
    
    printf("\nrun time (in seconds): %f on %d processors\n",t,p);
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    /* Print time */
    printf("\nrun time (in seconds): %f\n",t);

    return 0;
}
