/* Parallel Solutiuon to Steady-State heat Problem 
   Assumes that p divides n.
   Davis Knox
*/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <mpi.h>
//#define N 500
//#define EPSILON 0.01

int main(int argc, char *argv[]) {
	double diff;	/* Change in value */
	int i, j, k;	
	double mean;	/* Avg boundary Value */
	double **u; 	/* Old values */
	double *uStorage;
	double **w; 	/* New Values */
	double *wStorage;
	int N;
	double EPSILON;
	int print;
	int id;		// Process rank
	int p; 		// Number of processes
	int firstRow;   // index of first row (does not include ghost rows)
	int lastRow;    // index of last row (does not include ghost rows)
	int numRows;     // total number of rows, including ghost.
	MPI_Status stat;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	// Get stuff from command line
	if(argv[3] != NULL) {
		N = atoi(argv[1]);		// Grid size - int
		EPSILON = (double)atof(argv[2]);// EPSILON - double
		print = atoi(argv[3]);  	// Print? 1: 0;
		
	} else {
		printf("Incorrect Command Line Input... Terminating");
		exit(1);
	}
	
	// calculate number of rows per processor
	firstRow = id * N / p;
	lastRow = ((id + 1) * N) / p - 1;
	numRows = lastRow - firstRow + 1;
	
	if(!id || id == p - 1) {
		numRows += 1;
	} else {
		numRows += 2;
	}

	// Memory allocation
	uStorage = (double *) malloc(numRows * N * sizeof(double));
	u = (double **) malloc(numRows * sizeof(double*));
	
	// Do we have it?
	if(uStorage == NULL || u == NULL) {
		printf("Cannot Allocate Memory... Terminating");
		exit(1);
	}
	
	// set up array
	for(i = 0; i < numRows; i++) 
		u[i] = &uStorage[i * N];

	wStorage = (double *) malloc(numRows * N * sizeof(double));
	w = (double **) malloc(numRows * sizeof(double*));

	// Do we have it?
	if(wStorage == NULL || w == NULL){
		printf("Cannot Allocate Memory... Terminating");
		exit(1);
	}
	
	// set up array
	for(i = 0; i < numRows; i++) 
		w[i] = &wStorage[i * N];

	// start timer
	clock_t time = -clock();

	// Set boundary values and compute mean boundary value 
	mean = 300.0 / 4.0;

	// Init top boundary
	if(!id) {
		for(i = 0; i < N; i++) u[0][i] = 100.0;
	}

	// Init side boundaries
	for(i = 0; i < numRows; i++) {
		u[i][0] = 100.0;
		u[i][N-1] = 100.0;
	}

	// Init bottom boundary
	if(id == p - 1) {
		for(i = 1; i < N; i++) u[numRows-1][i] = 0.0;
	}

	// Init interior values
	for(i = 1; i < numRows-1; i++)
		for(j = 1; j < N-1; j++)
			u[i][j] = mean;

	// Init non end processor ends
	if(!(id == 0 || id == p - 1)) 
		for(i = 1; i < N-1; i++) {
			u[0][i] = mean;
			u[numRows-1][i] = mean;
		}

	
	double pDiff;
	// Compute steady-state solution
	for(;;){		
		diff = 0.0;
		for(i = 1; i < numRows-1; i++){
			for(j = 1; j < N-1; j++){
				w[i][j] = (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]) / 4.0;
				if(fabs(w[i][j] - u[i][j]) > diff){
					diff = fabs(w[i][j] - u[i][j]);
				}
			}
		}
		
		// Find largest diff to see if we continue
		MPI_Allreduce(&diff, &pDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); 
		

		// Stop if max diff is less than or equal to Epsilon
		if(pDiff <= EPSILON) break;

		for(i = 1; i < numRows-1; i++)
			for(j = 1; j < N-1; j++)
				u[i][j] = w[i][j];
		
		// Pass around ghost rows :: start with sending bottom rows
		if(id != p - 1) MPI_Send(u[numRows-2], N-1, MPI_DOUBLE, id+1, 4, MPI_COMM_WORLD);
		
		// Receive top rows
		if(id != 0) MPI_Recv(u[0], N-1, MPI_DOUBLE, id-1, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);

		// Send top rows
		if(id != 0) MPI_Send(u[1], N-1, MPI_DOUBLE, id-1, 4, MPI_COMM_WORLD);

		// Receive bottom rows
		if(id != p - 1) MPI_Recv(u[numRows-1], N-1, MPI_DOUBLE, id+1, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
			
	}
	
	// stop timer
	time += clock();

	// Print time
	double elapsed = (double)time / (double)CLOCKS_PER_SEC;
	if(!id) printf("Epsilon: %f\nN: %d\nProcessors: %d\nTime: %f seconds\n", EPSILON, N, p, elapsed);	


	if(print){
		// top boundary
		if(!id)for(i = 0; i < N; i++){
			printf("%6.2f ", u[0][j]);
		}
		putchar('\n');
		for(k = 0; k < p; k++){
			MPI_Barrier(MPI_COMM_WORLD);
			if(id == k){
				for(i = 1; i < numRows-1; i++) {
					for(j = 0; j < N; j++) {
						printf("%6.2f ", u[i][j]);
					}
					putchar ('\n');
				}
			}
		}
		// bottom boundary
		if(id == p - 1)for(i = 0; i < N; i++) {
			printf("%6.2f ", u[numRows-1][j]);
		}
		putchar('\n');
	}

	
	MPI_Finalize();

	return 0;
}
