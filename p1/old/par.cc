#include <math.h>
#include <cstdlib>
#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char* argv[])
{
	
	if (argc != 3) {
		cout << "Incorrect number of arguments. Should be: ssh.par [N] [EPSILON]" << endl;
		exit (1);
	}

	int id, p;
	long N;							/* N^2 blocks in the grid */
	double EPSILON;					/* minimum completion difference */
	
	int underEPSILON;				/* tracks if cycle was under EPSILON difference */
	double mean;					/* mean value to start in interrior blocks */
	double **oldv, **newv;			/* blocks held */
	double *ghost_n, *ghost_s;		/* north and south ghost rows */
	double elapsedTime;				/* elapsed time of the program */
	
	int low_row, high_row, local_rows;
	
	MPI::Init(argc, argv);
	
	id = MPI::COMM_WORLD.Get_rank();
	p = MPI::COMM_WORLD.Get_size();
	MPI_Barrier(MPI_COMM_WORLD);
	
	N = atoi(argv[1]);
	EPSILON = atof(argv[2]);
	
	//block-row decomposition...
	//		- process p gets [(n/p) x n] blocks (array slots)
	//			- p also gets 2 rows of ghost columns (except row 0 and size-1)
	//			- p communicates ghost columns with p-1 and p+1 (except for p=0,size-1
	
	
	// calc scatter/distribute values
	//		[integer division cover's floor operation]
	low_row = (long long) id*N /p;
	high_row = (long long) (id+1)*N /p - 1;
	local_rows = high_row - low_row + 1;

	// allocate memory for local grid and ghost rows
	//		- u, w will be n/p
	double *uSpace, *wSpace;
	uSpace = (double *) malloc (local_rows*N * sizeof(double));
	wSpace = (double *) malloc (local_rows*N * sizeof(double));
	oldv = (double **) malloc (local_rows * sizeof(double *));
	newv = (double **) malloc (local_rows * sizeof(double *));
	ghost_n = (double *) malloc (N * sizeof(double));	/* only first row doesn't have this */
	ghost_s = (double *) malloc (N * sizeof(double));	/* only last row doesn't have this  */
	
	if (uSpace==NULL || wSpace==NULL || oldv==NULL || newv==NULL || ghost_n==NULL || ghost_s==NULL)
	{
		cout << "Cannot allcate memory. Exiting." << endl;
		MPI::Finalize();
		exit (1);
	}
	
	for (int i=0; i<local_rows; i++) {
		oldv[i] = &uSpace[i*N];
		newv[i] = &wSpace[i*N];
	}
	
	/* Start Timing */
	elapsedTime = -MPI::Wtime();
	
	// initiate boarder values and mean values
	// if id==0, top row 100,
	//    id==p-1, bottom row 0.
	// all have side columns 100
	mean = (300.0 * N) / (4.0 * N);
	
	for (int i=0; i<local_rows; i++) {
		for (int j=0; j<N; j++) {
			if ((i==0 && low_row==0) || j==0 || j==N-1) {
				oldv[i][j] = newv[i][j] = 100.0;
			} else if (i==local_rows-1 && high_row==N-1) {
				oldv[i][j] = newv[i][j] = 0.0;
			} else {
				oldv[i][j] = mean;
			}
		}
	}
	
	// compute steady state distribution:
	int start_i = (low_row==0)? 1 : 0;
	int  end_i  = (high_row==N-1)? local_rows-2 : local_rows-1;
	underEPSILON = 0;
	
	while (!underEPSILON) {
		
		underEPSILON = 1;
		
		// send/refresh ghost rows
		if (id != p-1) {
			// send out south row to id+1
			MPI::COMM_WORLD.Send(oldv[local_rows-1], N, MPI::DOUBLE, id+1, 33);
		}
		if (id != 0)
		{
			// receive north ghost row from id-1
			MPI::COMM_WORLD.Recv(ghost_n, N, MPI::DOUBLE, id-1, MPI::ANY_TAG);
			
			// send out north row to id-1
			MPI::COMM_WORLD.Send(oldv[0], N, MPI::DOUBLE, id-1, 33);
		}
		if (id != p-1) {
			// receive south ghost row from id+1
			MPI::COMM_WORLD.Recv(ghost_s, N, MPI::DOUBLE, id+1, MPI::ANY_TAG);
		}
		
	
		// calculate next values for local grid
		for (int i=start_i; i <= end_i; i++)
		{
			for (int j=1; j < N-1; j++)
			{
				if (i == 0) // then we need to look at north ghost row
				{
					newv[i][j] = ( ghost_n[j]   + oldv[i][j-1] +
								   oldv[i+1][j] + oldv[i][j+1] ) / 4.0;
				}
				else if (i == local_rows-1) // then we need to look at south ghost row
				{
					newv[i][j] = ( oldv[i-1][j] + oldv[i][j-1] +
								   ghost_s[j]   + oldv[i][j+1] ) / 4.0;
				}
				else
				{
					newv[i][j] = ( oldv[i-1][j] + oldv[i][j-1] +
								   oldv[i+1][j] + oldv[i][j+1] ) / 4.0;
				}
				
				// if the current block changed under threashold
				if (fabs(newv[i][j] - oldv[i][j]) > EPSILON)
					underEPSILON = 0;
			}
		}
		
		//need all-to-all reduction of epsilon bool
		MPI_Allreduce(&underEPSILON, &underEPSILON, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
		
		// swap new and old array pointers
		double **tmp = newv;
		newv = oldv;
		oldv = tmp;
	}
	
	
	// reduction to root here?
	
	
	// stop timing
	elapsedTime += MPI::Wtime();
	
	if (!id) cout << endl << "processors = " << p << " time: " << elapsedTime << endl << endl;
	
	// print out values to ensure correctness
	// let each processor print out what it has, one processor at a time
	for (int ident=0; ident<p; ident++) {
		MPI_Barrier(MPI_COMM_WORLD);
		
		if (id == ident)
		{
			for (int i=0; i<local_rows; i++)
			{
				cout << "R[" << low_row+i << "]: \t";
				for (int j=0; j<N; j++)
					printf("%6.2f ", oldv[i][j]);
				cout << endl;
			}
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI::Finalize();
	return 0;

}
