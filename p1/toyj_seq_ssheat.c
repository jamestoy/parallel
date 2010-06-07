/* Sequential Solution to Steady-State Heat Problem */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

int
main (int argc, char *argv[])
{
    int n;
    double eps;
    if(argc != 3) {
        printf("usage: ssheat <grid size> <epsilon>\n");
        return -1;
    } else {
        n = atoi(argv[1]);
        eps = atof(argv[2]);
    }
    double diff;    /* Change in value */
    int i,j,iter;
    double mean;    /* Average boundary value */
    double **u;
    double **w;

    assert((u = malloc(n * sizeof(double *))) != NULL);
    assert((w = malloc(n * sizeof(double *))) != NULL);
    for(i = 0; i < n; i++) {
        assert((u[i] = malloc(n * sizeof(double))) != NULL);
        assert((w[i] = malloc(n * sizeof(double))) != NULL);
    }
    
    clock_t start, stop;
    double t = 0.0;
    assert((start = clock()) != -1);

    /* Set boundary values and compute mean boundary value */
    mean = 0.0;
    for (i = 0; i < n; i++) {
        u[i][0] = u[i][n-1] = u[0][i] = 100.0;
        u[n-1][i] = 0.0;
        mean += u[i][0] + u[i][n-1] + u[0][i] + u[n-1][i];    
        
        //w[i][0] = w[i][n-1] = w[0][i] = 100.0;
        //w[n-1][i] = 0.0;
    }
    mean /= (4.0 * n);

    /* Initialize interior values */
    for(i = 1; i < n-1; i++)
        for(j = 1; j < n-1; j++) u[i][j] = mean;

    /* Compute the steady-state solution */
    iter = 0;
    for(;;) {
        diff = 0.0;
        for(i = 1; i < n-1; i++)
            for(j = 1; j < n-1; j++) {
                w[i][j] = (u[i-1][j] + u[i+1][j] + 
                           u[i][j-1] + u[i][j+1]) / 4.0;
                if (fabs(w[i][j] - u[i][j]) > diff)//XXX: changed to j+1....
                    diff = fabs(w[i][j] - u[i][j]);
            }
        iter++;
        if (diff <= eps) break;

        for(i = 1; i < n-1; i++)
            for(j = 1; j < n-1; j++) u[i][j] = w[i][j];
    }
    /* Print number of iteration */
    printf("number of iterations: %d\n\n",iter);
    
    /* Derive run-time */
    stop = clock();
    t = (double)(stop-start)/CLOCKS_PER_SEC;

    /* Print solution */
    for (i = 0; i < n; i++) {
        for(j = 0; j < n; j++)
            printf("%6.2f ", w[i][j]);
        putchar('\n');
    }

    /* Print time */
    printf("\nrun time (in seconds): %f\n",t);
}
