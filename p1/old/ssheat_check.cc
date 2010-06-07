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
	int p_side;
	long N;							/* N^2 blocks in the grid */
	double EPSILON;					/* minimum completion difference */
	
	int underEPSILON;				/* tracks if cycle was under EPSILON difference */
	double mean;					/* mean value to start in interrior blocks */
	double **oldv, **newv;			/* blocks held */
	double *ghost_n, *ghost_s,
		   *ghost_e, *ghost_w;		/* north and south ghost rows */
	double *tmp;
	double elapsedTime;				/* elapsed time of the program */
	
	int low_row, high_row, low_col, high_col, local_rows, local_cols;
	
	MPI::Init(argc, argv);
	
	id = MPI::COMM_WORLD.Get_rank();
	p = MPI::COMM_WORLD.Get_size();
	p_side = (int)sqrt(p);
	MPI_Barrier(MPI_COMM_WORLD);
	
	N = atoi(argv[1]);
	EPSILON = atof(argv[2]);
	
	if ( sqrt(p) - (int)(sqrt(p)) != 0.00)
	{
		cout << "ERROR :: Number of processors not fit for a square grid." << endl;
		MPI::Finalize();
		exit (1);
	}
	
	// calc scatter/distribute values
	//		[integer division cover's floor operation]
	low_row = (int)floor((id / p_side) * (N / (double)p_side));
	high_row = (int)floor(((id / p_side) + 1) * (N / (double)p_side)) - 1;
	local_rows = high_row - low_row + 1;
	
	low_col = (int)floor((id % p_side) * (N / (double)p_side));
	high_col = (int)floor(((id % p_side) + 1) * (N / (double)p_side)) - 1;
	local_cols = high_col - low_col + 1;
	
	//TEST
//	cout << "ID("<<id<<") :: Row_f="<<low_row<<", Row_l="<<high_row<<", Col_f="<<low_col<<", Col_l="<<high_col<<endl;

	// allocate memory for local grid and ghost rows
	//		- u, w will be n/p
	double *uSpace, *wSpace;
	uSpace = (double *) malloc (local_rows*local_cols * sizeof(double));
	wSpace = (double *) malloc (local_rows*local_cols * sizeof(double));
	oldv = (double **) malloc (local_rows * sizeof(double *));
	newv = (double **) malloc (local_rows * sizeof(double *));
	ghost_n = (double *) malloc (local_cols * sizeof(double));	/* only first row doesn't use this */
	ghost_s = (double *) malloc (local_cols * sizeof(double));	/* only last  row doesn't use this */
	ghost_e = (double *) malloc (local_rows * sizeof(double));	/* only the rightmost column doesn't use this */
	ghost_w = (double *) malloc	(local_rows * sizeof(double));	/* only the leftmost  column doesn't use this */
	tmp = (double *) malloc (local_rows * sizeof(double));
	
	if (   uSpace==NULL || wSpace==NULL
		|| oldv==NULL || newv==NULL
		|| ghost_n==NULL || ghost_s==NULL
		|| ghost_e==NULL || ghost_w==NULL || tmp==NULL)
	{
		cout << "Cannot allcate memory. Exiting." << endl;
		MPI::Finalize();
		exit (1);
	}
	
	for (int i=0; i<local_rows; i++) {
		oldv[i] = &uSpace[i*local_cols];
		newv[i] = &wSpace[i*local_cols];
	}
	
	/* Start Timing */
	elapsedTime = -MPI::Wtime();
	
	// initiate boarder values and mean values
	// if id==0, top row 100,
	//    id==p-1, bottom row 0.
	// all have side columns 100
	mean = (300.0 * (N-2)) / (4.0 * (N-2)); //average will be 75.00 every time, so...
	
	for (int i=0; i<local_rows; i++) {
		for (int j=0; j<local_cols; j++) {
			if (( i == 0  && low_row == 0 )
			 || ( j == 0  && low_col == 0 )
			 || ( j == local_cols-1 && high_col == N-1 )) {
				oldv[i][j] = newv[i][j] = 100.0;
			} else if (i==local_rows-1 && high_row==N-1) {
				oldv[i][j] = newv[i][j] = 0.0;
			} else {
				oldv[i][j] = mean;
			}
		}
	}
	
	// compute steady state distribution:
	underEPSILON = 0;
	int start_i = (low_row==0)? 1 : 0;
	int  end_i  = (high_row==N-1)? local_rows-2 : local_rows-1;
	int start_j = (low_col==0)? 1 : 0;
	int  end_j  = (high_col==N-1)? local_cols-2 : local_cols-1;
	int i;
	double n, s, e, w; /* hold values of north, south, eart and west blocks */
	
	while (!underEPSILON) {
		
		underEPSILON = 1;
		
		
		// send/receive north/south ghost rows
		if (id/p_side != p_side-1)	// not the bottom row
		{
			MPI::COMM_WORLD.Send(oldv[local_rows-1], local_cols, MPI::DOUBLE, id+p_side, 44);	//SOUTH side SEND down
		}
		if (id/p_side != 0)			// not the top row
		{
			MPI::COMM_WORLD.Recv(ghost_n, local_cols, MPI::DOUBLE, id-p_side, MPI::ANY_TAG);	//NORTH side RECEIVE up
			MPI::COMM_WORLD.Send(oldv[0], local_cols, MPI::DOUBLE, id-p_side, 44);				//NORTH side SEND up
		}
		if (id/p_side != p_side-1)	// not the bottom row
		{
			MPI::COMM_WORLD.Recv(ghost_s, local_cols, MPI::DOUBLE, id+p_side, MPI::ANY_TAG);	//SOUTH side RECEIVE down
		}
		
		// send/receive east/west ghost rows
		if (id%p_side != p_side-1)
		{
			//use ghost_w to build east side to send, SEND to RIGHT
			for ( i=0; i<local_rows; i++)
				tmp[i] = oldv[i][local_cols-1];
			MPI::COMM_WORLD.Send(tmp, local_rows, MPI::DOUBLE, id+1, 44);
		}
		if (id%p_side != 0)
		{
			//RECEIVE from LEFT, insert to ghost_w
			MPI::COMM_WORLD.Recv(ghost_w, local_rows, MPI::DOUBLE, id-1, MPI::ANY_TAG);
			
			//use ghost_e to build west side to send, SEND to LEFT
			for ( i=0; i<local_rows; i++)
				tmp[i] = oldv[i][0];
			MPI::COMM_WORLD.Send(tmp, local_rows, MPI::DOUBLE, id-1, 44);
		}
		if (id%p_side != p_side-1)
		{
			//RECEIVE from RIGHT, insert to ghost_e
			MPI::COMM_WORLD.Recv(ghost_e, local_rows, MPI::DOUBLE, id+1, MPI::ANY_TAG);
		}
		
		
		// calculate next values for local grid
		for (int i=start_i; i <= end_i; i++)
		{
			for (int j=start_j; j <= end_j; j++)
			{
				//north block
				if (i==0) n = ghost_n[j];
				else n = oldv[i-1][j];
				
				//south block
				if (i==local_rows-1) s = ghost_s[j];
				else s = oldv[i+1][j];
				
				//east block
				if (j==local_cols-1) e = ghost_e[i];
				else e = oldv[i][j+1];
				
				//west block
				if (j==0) w = ghost_w[i];
				else w = oldv[i][j-1];
				
				newv[i][j] = (n+s+e+w) / 4.0;
				
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
	
	if (!id) cout << endl << "Time Elapsed over " << p << " processors = " << elapsedTime << endl << endl;
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI::Finalize();
	return 0;

}
