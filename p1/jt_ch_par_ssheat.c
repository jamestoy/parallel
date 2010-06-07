/* Parallel (checkerboard decomposition) Solution to Steady-State Heat Problem */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>

#define TAG 3
#define INIT_MEAN 75.0
#define NIL 0

/* PROTOTYPES */
void print_solution(double **u,int id,int p,int localc,int localr,int n);
int compute_solution(int _p, int id, double eps, 
					 int startr, int endr, int startc, int endc,
					 int localc, int localr, double *_north, double *_east, 
					 double *_south, double *_west, double **u, double **w, double *row_hold);

int compute_solution(int _p, int id, double eps, 
						int startr, int endr, int startc, int endc,
						int localc, int localr, double *_north, double *_east, 
						double *_south, double *_west, double **u, double **w, double *row_hold)
{
	
	double diff,
			tdiff,
			global_diff;     /* Change in value */
	
	double north,
			south,
			east,
			west;				/* Holds specific block value */
	
	int i,j,iter;
	MPI_Status stat;
	
	/* Compute the steady-state */
	iter = 0;
    for (;;) {
		/* east and west */
		if (id/_p != _p-1)//3
			MPI_Recv(_south, localc, MPI_DOUBLE, id+_p, TAG, MPI_COMM_WORLD, &stat);
		if (id%_p != _p-1)//6
			MPI_Recv(_east, localr, MPI_DOUBLE, id+1, TAG, MPI_COMM_WORLD, &stat);
		if (id%_p != _p-1){//4
			for (i = 0; i < localr; i++)
				row_hold[i] = u[i][localc--];
			MPI_Send(row_hold, localr, MPI_DOUBLE, id+1, TAG,MPI_COMM_WORLD);
		}
		if (id%_p != 0) {//5
			MPI_Recv(_west, localr, MPI_DOUBLE, id-1, TAG, MPI_COMM_WORLD, &stat);
			for (i = 0; i < localr; i++)
				row_hold[i] = u[i][0];
			MPI_Send(row_hold, localr, MPI_DOUBLE, id-1, TAG, MPI_COMM_WORLD);
		}
		/* north and south */
		if (id/_p != _p-1)//1 
			MPI_Send(u[localr-1], localc, MPI_DOUBLE, id+_p, TAG, MPI_COMM_WORLD);
		if (id/_p != 0) { //2
			MPI_Recv(_north, localc, MPI_DOUBLE, id-_p, TAG, MPI_COMM_WORLD, &stat);
			MPI_Send(u[0], localc, MPI_DOUBLE, id-_p, TAG, MPI_COMM_WORLD);
		}
		
		/* Set diff back to 0 and generate values for each grid */
		diff = 0.0;
        for(i = startr; i <= endr; i++)
            for(j = startc; j <= endc; j++) {
				/* Never Eat Shredded Wheat */
				north = (i == 0) ? _north[j] : u[i-1][j];
				east = (j == localc-1) ? _east[i] : u[i][j+1];
				south = (i == localr-1) ? _south[j] : u[i+1][j];
				west = (j == 0) ? _west[i] : u[i][j-1];
				w[i][j] = (north + east + south + west) / 4;
				/* Delimeter */
				if (fabs(w[i][j] - u[i][j]) > diff)
					diff = fabs(w[i][j] - u[i][j]);
            }
		double **oldptr = w;
		w = u;
		u = oldptr;
		
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
   
	int i,j,its,           // iterators
		 startr,endr,
        startc,endc,
        lr,hr,localr,       // rows
        lc,hc,localc,       // columns
        id,p,_p;            // constants
   
	double mean;            /* Average boundary value */
    
	double **u;				/* Old values */
	double **w;				/* New values */
	double *_north,			/* Ghost rows */
			*_south,
			*_east,
			*_west,
			*row_hold;			/* For building ghost rows */
	
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
	
	/* Make sure we have square proc grid */
	assert((sqrt(p)-sqrt(p)) == 0);
	
	/* Delimits side of a proc */
	_p = sqrt(p);
	
    /* Scatter from class 
     * low row, high row, local row
     *
     * XXX: (jft 20100521) added columns for checkerboard
     */
    lr = id * n / p; lc = (id % _p) * (n / _p);
    hr = (id + 1) * n /p - 1; hc = (id % _p + 1) * (n / _p) - 1;
    localr = hr - lr + 1; localc = hc - lc + 1;

    /* Allocation step */
    assert((u = malloc(localr * localc * sizeof(*u))) != NULL);
    assert((w = malloc(localr * localc * sizeof(*w))) != NULL);
    for(i = 0; i < localr; i++) {
        //XXX: (jft) 20100520 not sure if i can do 
		 //sizeof(**u) & sizeof(**w) seems to be working
        assert((u[i] = malloc(localc * sizeof(**u))) != NULL);
        assert((w[i] = malloc(localc * sizeof(**w))) != NULL);
    }
	assert((_north = malloc(localr * sizeof(*_north))) != NULL);
	assert((_south = malloc(localr * sizeof(*_south))) != NULL);
	assert((_east = malloc(localc * sizeof(*_east))) != NULL);
	assert((_west = malloc(localc * sizeof(*_west))) != NULL);
	/* could be localc or localr */
	assert((row_hold = malloc(localr * sizeof(*row_hold))) != NULL);
    
   /* Prepare to be timed */ 
	double t = MPI_Wtime();
    
    /* Set boundary values and compute mean boundary value
     * XXX: (jft) 20100521 - need to set limits based on rows and cols
     *                       otherwise everything will get clobbered
     */

    for(i = 0; i < localr; i++)
        for(j = 0; j < localc; j++) {
            if(i==localr-1 && hr==n-1) //bottom 
                w[i][j] = u[i][j] = 0.0;
            else if((j==localc-1 && hc==n-1) || (i==0 && lr==0) || (j==0 && lc==0))
                w[i][j] = u[i][j] = 100.0;
            else
                u[i][j] = INIT_MEAN;
        }
    
	/* Setup checkerboard boundaries based on id / p */
    if(lr == 0) startr = 1;             // start out 1
    else startr = 0;                    // go all the way
    if(hr == n -1) endr = localr - 2;   // two out
    else endr = localr - 1;             // go to end

    if(lc == 0) startc = 1;
    else startc = 0;
    if(hc == n-1) endc = localc - 2;
    else endc = localc - 1;
	
	/* Compute Solution */
	its = compute_solution(_p,id,eps,startr,endr,startc,endc,
						   localc,localr,_north,_east,_south,_west,u,w,row_hold);
	
    /* Print number of iteration */
    printf("number of iterations: %d\n",its);
    
    /* Derive run-time */
    t += MPI_Wtime();
    
	/* Print time */
	printf("run time (in seconds): %f on %d processors\n",t,p);
	
	/* Deallocation */
	free(u); free(w); free(_north); free(_south);
    
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	
	return 0;
}
