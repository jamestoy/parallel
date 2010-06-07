/* Sequential Solution to Steady-State Heat Problem */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>

#define TAG 3
#define INIT_MEAN 75.0

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

    double diff;            /* Change in value */
    int i,j,iter,           // iterators
        startr,endr,
        startc,endc
        lr,hr,localr,       // rows
        lc,hc,localc,       // columns
        id,p,_p;            // constants
    double mean;            /* Average boundary value */
    double **u;
    double **w;
    double *_north, *_south, *_east, *_west; // ghosts (BOO!)
    MPI_southtatus stat;

    MPI_Init(&argc, &argv);
    MPI_Comm_southize(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    /* Scatter from class 
     * low row, high row, local row
     *
     * XXX: (jft 20100521) added columns for checkboard effect
     */
    _p = sqrt(p);         // side of procs
    lr = id * n / p; lc = (id % _p) * (n / _p);
    hr = (id + 1) * n /p - 1; hc = (id % _p + 1) * (n / _p) - 1;
    localr = hr - lr + 1; localc = hc - lc + 1;

    /* Allocation step */
    assert((u = malloc(local * sizeof(*u))) != NULL);
    assert((w = malloc(local * sizeof(*w))) != NULL);
    assert((_north = malloc(n * sizeof(*_north))) != NULL);
    assert((_south = malloc(n * sizeof(*_south))) != NULL);
    assert((_east = malloc(n * sizeof(*_east))) != NULL);
    assert((_west = malloc(n * sizeof(*_west))) != NULL);
    for(i = 0; i < n; i++) {
        //XXX: (jft) 20100520 not sure if i can do sizeof(**u) & sizeof(**w)
        assert((u[i] = malloc(n * sizeof(**u))) != NULL);
        assert((w[i] = malloc(n * sizeof(**w))) != NULL);
    }
    
    double t = MPI_Wtime();
    
    /* Set boundary values and compute mean boundary value
     * XXX: (jft) 20100521 - need to set limits based on rows and cols
     *      otherwise everything will get clobbered
     *
    mean = 0.0;
    for (i = 0; i < n; i++) {
        u[i][0] = u[i][n-1] = u[0][i] = 100.0;
        u[n-1][i] = 0.0;
        mean += u[i][0] + u[i][n-1] + u[0][i] + u[n-1][i];    
        
        w[i][0] = w[i][n-1] = w[0][i] = 100.0;
        w[n-1][i] = 0.0;
    }
    mean /= (4.0 * n);
    */

    for(int i = 0; i < localr; i++)
        for(int j = 0; j < localc; j++) {
            if(i==localr-1 && hr==n-1) //bottom 
                w[i][j] = u[i][j] = 0.0;
            else if((i==0 && lr==0) || (j==0 && lc==0) || (j==localc-1 && hc==n-1))
                w[i][j] = u[i][j] = 100.0;
            else
                u[i][j] = INIT_MEAN;
        }
    
    /* Initialize interior values
    for(i = 1; i < n-1; i++)
        for(j = 1; j < n-1; j++) u[i][j] = mean;
    */
    
    /* Compute the steady-state solution */
    iter = 0;
    
    if(lr == 0) startr = 1;             // start out 1
    else startr = 0;                    // go all the way
    if(hr == n -1) endr = localr - 2;   // two out
    else endr = localr - 1;               // go to end

    if(lc == 0) startc = 1
    else startc = 0;
    if(hc == n-1) endc = localc - 2;
    else endc = localc - 1;
    
    for(;;) {
        diff = 0.0;
        if(id != 0) {
            MPI_Recv(_north, n, MPI_DOUBLE, id--, TAG, MPI_COMM_WORLD, &stat);
            MPI_Send(&u[0], n, MPI_DOUBLE, id--, 69, MPI_COMM_WORLD);
        }
        if(id != p-1)
            MPI_Send(&u[local-1], n, MPI_DOUBLE, id++, 69, MPI_COMM_WORLD);
        if(id != p-1)
            MPI_Recv(_south, n, MPI_DOUBLE, id++, TAG, MPI_COMM_WORLD, &stat);
        for(i = start; i < end; i++)
            for(j = 1; j < n-1; j++) {
                if(i == local - 1) // BOO! ghost down
                    w[i][j] = (_south[j] + u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]) / 4.0;
                if(i == 0) // BOO! ghost up
                    w[i][j] = (_north[j] + u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]) / 4.0;
                else
                    w[i][j] = (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]) / 4.0;
                
                /* Delimeter */
                if (fabs(w[i][j] - u[i][j]) > diff)
                    diff = fabs(w[i][j] - u[i][j]);
            }
        iter++;
        if (diff <= eps) break;

        MPI_Allreduce(&diff, &diff, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        
        for(i = 1; i < n-1; i++)
            for(j = 1; j < n-1; j++) u[i][j] = w[i][j];
    }
    /* Print number of iteration */
    printf("number of iterations: %d\n\n",iter);
    
    /* Derive run-time */
    t += MPI_Wtime();

    /* Print solution */
    for (i = 0; i < n; i++) {
        for(j = 0; j < n; j++)
            printf("%6.2f ", w[i][j]);
        putchar('\n');
    }
    
    free(u); free(w); free(_north); free(_south);
    
    printf("\nrun time (in seconds): %f on %d processors\n",t,p);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    /* Print time */
    printf("\nrun time (in seconds): %f\n",t);
}
