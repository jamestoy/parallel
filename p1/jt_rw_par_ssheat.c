/* Parallel (row-wise decomposition) Solution to Steady-State Heat Problem.
 * this implementation is similar to the one found in chapter 18 of the 
 * Quinn text.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>

#define TAG 0
#define NIL 0

/* PROTOTYPES */
void print_solution(double **u, int id, int p, int local, int n);
int compute_solution(int p, int id, int n, double eps,
					    int start, int end, int local, 
					    double *_n, double *_s, double **u, double **w);

/* Print solution */
void print_solution(double **u, int id, int p, int local, int n)
{
	int i,j;
	int iter_id;
	
	for (iter_id = 0; iter_id < 0; iter_id++) {
		MPI_Barrier(MPI_COMM_WORLD);
		for (i = 0; i < local; i++) {
			for(j = 0; j < n; j++)
				printf("%6.2f ", u[i][j]);
			putchar('\n');
		}
	}
}

int compute_solution(int p, int id, int n, double eps,
					    int start, int end, int local, 
					    double *_n, double *_s, double **u, double **w)
{
	double diff,
	tdiff,
	global_diff;    /* Change in value */
	
	int i,j,iter;
	MPI_Status stat;
	
	/* Compute the steady-state solution */
    iter = 0;
    for (;;) {        
        if (id != p-1)
            MPI_Send(w[local-1], n, MPI_DOUBLE, id+1, TAG, MPI_COMM_WORLD);
        if (id != 0) {
            MPI_Recv(_n, n, MPI_DOUBLE, id-1, TAG, MPI_COMM_WORLD, &stat);
            MPI_Send(w[0], n, MPI_DOUBLE, id-1, TAG, MPI_COMM_WORLD);
        }
        if (id != p-1)
            MPI_Recv(_s, n, MPI_DOUBLE, id+1, TAG, MPI_COMM_WORLD, &stat);
		diff = 0.0;
//#pragma omp parallel private(i,j,tdiff)
//		{
//			tdiff = 0.0;
//#pragma omp for
			for(i = start; i <= end; i++)    //XXX:(jft) s/</<=/g
				for(j = 1; j < n-1; j++) {
					if (i == 0) {
						w[i][j] = (_n[j] + 
								    u[i][j+1] + 
								    u[i][j-1] + 
								    u[i+1][j]) / 4.0;
					
					} else if (i == local - 1) {
						w[i][j] = (_s[j] + 
								    u[i][j+1] + 
								    u[i][j-1] + 
								    u[i-1][j]) / 4.0;
					} else {                    // regular computation
						w[i][j] = (u[i][j+1] + 
								    u[i][j-1] + 
								    u[i+1][j] + 
								    u[i-1][j]) / 4.0;
					}
					
					/* Delimeter */
					if (fabs(w[i][j] - u[i][j]) > diff) //XXX: (jft) s/diff/tdiff
						diff = fabs(w[i][j] - u[i][j]);
				}
//#pragma omp for nowait
			for(i = 1; i < n-1; i++)
				for(j = 1; j < n-1; j++)
					u[i][j] = w[i][j];
//#pragma omp critical
//			if (tdiff > diff) diff = tdiff;
//		}
        MPI_Allreduce(&diff, &global_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (global_diff <= eps) break;
        iter++;
    }
	return iter;
}

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
    
    int i,j,its,
        lr,hr,local,
        id,p,
        start,end;

    double mean;                 /* Average boundary value */
    
    double **u; /* Old values */
    double **w; /* New values */
    double *_n, /* North ghost */
           *_s; /* South ghost */

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
    assert((u = malloc(local * sizeof(double *))) != NULL);
    assert((w = malloc(local * sizeof(double *))) != NULL);
    for(i = 0; i < local; i++) { //XXX:(jft) s/n/local
        assert((u[i] = malloc(n * sizeof(double))) != NULL);
        assert((w[i] = malloc(n * sizeof(double))) != NULL);
    }
    /* Ghost Rows */
    assert((_n = malloc(n * sizeof(double))) != NULL);
    assert((_s = malloc(n * sizeof(double))) != NULL);

    double t = MPI_Wtime();
    
    /* Set boundary values and compute mean boundary value */
    mean = 75.0;
    for (i = 0; i < local; i++) {
		for (j = 0; j < n; j++) {
			if (j == n-1 || j == 0 || (i == NIL && lr == NIL)) { // LEFT / TOP / RIGHT
                u[i][j] = w[i][j] = 100.0;
            } else if (i == local-1 && hr == n-1) {				 // BOTTOM
                u[i][j] = w[i][j] = 0.0;
            } else {
                u[i][j] = mean;									 // ALL INTERIOR VALUES TO 75.0
            }
        }
    }

    /* Set iterator boundaries */
    if(lr == 0)
        start = 1;          
    else 
        start = 0;          
    if(hr == n-1) 
        end = local-2;
    else
        end = local-1;
	
	/* End of initialization */

	its = compute_solution(p,id,n,eps,
							  start,end,local,
							  _n,_s,u,w);
	
	print_solution(u,id,p,local,n);
	
    /* Print number of iterations */
    printf("number of iterations: %d\n\n",its);
    
    /* Derive run-time */
    t += MPI_Wtime();
    
    if (!id) printf("\nrun time (in seconds): %f on %d processors\n",t,p);
    
	/* Deallocation */
	free(u); free(w); free(_n); free(_s);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
