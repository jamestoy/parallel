/*
*Matrix-vector multiplication, Version 3
*
*This program multiplies a matrix and a vector.
*The matrix and vector are input from files.
*The answer is printed to standard output.
*Does only work for square p (number of processors).
*
*Data distribution of matrix: checkerboard block-wise
*Data distribution of vector: block
*
*
Programmed by Annika Biermann
*
*Last modification: 24 November 2009
*/

#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include "./MyMPI.h"
#include "./MyMPI.c"
/*Change these two definitions when the matrix and vector element types change */
typedef double dtype;
#define mpitype MPI_DOUBLE
int main (int argc, char *argv[]) {
dtype **a;	/* First factor, a matrix */
dtype *b;	/* Second factor, a vector */
dtype *v;	/* Local vector block */
dtype *c_part_out;	/* Partial sums, sent */
dtype *storage;	/* Matrix elements stored here */
int p;	/* Number of processes */
int size[2];	/* Vector with size of each grid dimension*/
int periodic[2];	/* Message wraparound flags*/
int id; 	/* Rank of process in virtual grid */
int i, j;	/* Loop indice */
int grid_coords[2];	/* Location of process in grid */
int m; /* Rows in the matrix */
int n;	/* Columns in the matrix */
int	nprime;	/* # of elements in the local initial vector block */
int nfinal;	/* # of elements in the local final vector block */
int dims;	/* # of dimensions of the grid */
int	local_cols; /* Local # of columns */
int	local_rows; /* Local # of rows */
double max_seconds;
double seconds; /* Elapsed time for matrix-vector multiply */
MPI_Comm grid_comm;	/* 2-D process grid, Cartesian topology	communicator */
MPI_Comm row_comm;	/* Processes in same row */
MPI_Comm column_comm;	/* Processes in same column */
MPI_Status status;	/* Result of receive */

/* Initialization */
MPI_Init (&argc, &argv);
MPI_Comm_rank (MPI_COMM_WORLD, &id);
MPI_Comm_size (MPI_COMM_WORLD, &p);
dims = 2;
nprime = 0;

/* Create virtual grid */
size[0] = size[1] = 0;
MPI_Dims_create(p, dims, size);
periodic[0] = periodic[1] = 0;
MPI_Cart_create(MPI_COMM_WORLD, dims, size, periodic, 1, &grid_comm);
if(size[0] == size[1]) {
	/* Read and print checkerboard matrix */
	read_checkerboard_matrix (argv[1], (void ***) &a, (void **) &storage, mpitype, &m, &n, grid_comm);
	print_checkerboard_matrix ((void **) a, mpitype, m, n, grid_comm);
	/* Read and print vector */
	MPI_Cart_coords (grid_comm, id, dims, grid_coords);
	MPI_Comm_split (grid_comm, grid_coords[1], grid_coords[0], &column_comm);
	if(grid_coords[1] == 0) { // only processes of the first column
		read_block_vector (argv[2], (void **) &b, mpitype, &nprime, column_comm);
		print_block_vector ((void *) b, mpitype, nprime,column_comm);
	}
	
	/* Redistribute vector */
	redistribute_vector((void *) b, (void **) &v, mpitype, n, nprime, (int *) &nfinal, size[0], size[1], grid_coords[0], grid_coords[1], grid_comm);
	/* Matrix-vector multiply */
	/* Each process multiplies its submatrix of 'a' and block of vector 'b', resulting in a partial block sum of product 'c'. */
	local_rows = BLOCK_SIZE(grid_coords[0], size[0], m);
	local_cols = BLOCK_SIZE(grid_coords[1], size[1], n);
	c_part_out = (dtype *) my_malloc (id, local_rows *sizeof(dtype));
	MPI_Barrier (MPI_COMM_WORLD);
	seconds = -MPI_Wtime();
	for (i = 0; i < local_rows; i++) {
		c_part_out[i] = 0.0;
		for (j = 0; j < local_cols; j++)
			c_part_out[i] += a[i][j] * v[j];
	}
	/* Reduce vectors across rows */
	/* Split grid_communicator into row_communicators */
	MPI_Comm_split (grid_comm, grid_coords[0], grid_coords[1], &row_comm);
	dtype c_part_in[size[1]-1][local_rows];

	/* If process not in first grid column */
	if (grid_coords[1] != 0) {
		/* Send block of c to first process in row */
		MPI_Send (c_part_out, local_rows, mpitype, 0, DATA_MSG, row_comm);
	}
	/* if process in first grid column */
	else {
		// for each column except the first one
		for (i = 1; i < size[1]; i++) {
			/* Receive blocks of c from other processes and do summation */
			MPI_Recv (c_part_in[i-1], local_rows, mpitype, i, DATA_MSG, row_comm, &status);
		}
		for (i = 0; i < local_rows; i++) {
			for(j=0;j<size[1]-1;j++){
				c_part_out[i] += c_part_in[j][i];
			}
		}
	}
	MPI_Barrier (MPI_COMM_WORLD);
	seconds += MPI_Wtime();
	MPI_Allreduce (&seconds, &max_seconds, 1, mpitype, MPI_MAX, MPI_COMM_WORLD);
	if (!id) {
		printf ("MV3) N = %d, Processes = %d, Time = %12.6f sec,", n, p, max_seconds);
		printf ("Mflop = %6.2f\n", 2*n*n/(1000000.0*max_seconds));
	}
	if(grid_coords[1] == 0)
		print_block_vector ((void *) c_part_out, mpitype, n, column_comm);
}
MPI_Finalize();
return 0;
}

/*
* This function is used to redistribute a vector block from the first process grid column processes
* to all processes within a communicator.
*/
void redistribute_vector(
void *b, 			/* IN - Address of vector */
void **v,			/* OUT - Subvector */
MPI_Datatype dtype,	/* IN - Vector element type */
int n,				/* IN - Size of (whole) input vector */
int nprime,			/* IN - Elements in initial block vector */
int *nfinal,		/* OUT - Elements in final block vector */
int grid_rows,		/* IN - Number of rows in the virtual process grid */
int grid_cols,		/* IN - Number of columns in the virtual process grid */
int grid_row_coord,	/* IN - row index of process in virtual process grid */
int grid_col_coord,	/* IN - column index of process in virtual process grid */
MPI_Comm grid_comm)	/* IN - Communicator */
{
int datum_size;		/* Bytes per vector element */
MPI_Status status;	/* Result of receive */
int id; 			/* Process rank */
int p;				/* Number of processes */
int dest_coords[2];	/* Coordinates of receiving process */
int	dest_id;		/* Rank of receiving process */
int sent_coords[2];	/* Coordinates of sending process */
int sent_id;		/* Rank of sending process */
MPI_Comm column_comm;/* Virtual grid column communicator */

MPI_Comm_size (grid_comm, &p);
MPI_Comm_rank (grid_comm, &id);
datum_size = get_size (dtype);
*nfinal = BLOCK_SIZE(grid_row_coord, grid_cols, n);

/* Allocate memory for final local vector block */
*v = my_malloc (id, (*nfinal) * datum_size);
/* if p is a square number (grid_rows == grid_cols) */
if(grid_rows == grid_cols) {
	/* Send/receive blocks of vector */
	/* if process of first grid column */
	if(grid_col_coord == 0) {
		/* if process 0 */
		if (grid_row_coord == 0) {
		/* copy vector block */
			memcpy(*v, b, (*nfinal) * datum_size);
		}
		else {
			/* Set coordinates */
			dest_coords[0] = grid_col_coord;
			dest_coords[1] = grid_row_coord;
			/* Get rank of destinating process */
			MPI_Cart_rank (grid_comm, dest_coords, &dest_id);
			/* Send vector block */
			MPI_Send (b, *nfinal, dtype, dest_id, DATA_MSG, grid_comm);
		}
	}
	/* if process of first grid row (except process 0) */
	else if (grid_row_coord == 0) {
		/* Set coordinates */
		sent_coords[0] = grid_col_coord;
		sent_coords[1] = grid_row_coord;
		/* Get rank of sending process */
		MPI_Cart_rank (grid_comm, sent_coords, &sent_id);
		/* Receive vector block */
		MPI_Recv (*v, *nfinal, dtype, sent_id, DATA_MSG, grid_comm, &status);
	}
	/* Broadcast blocks from first row per column */
	/* Split grid communicator into column communicator */
	MPI_Comm_split (grid_comm, grid_col_coord, grid_row_coord, &column_comm);
	/* Broadcast vector block per column */
	MPI_Bcast (*v, *nfinal, dtype, 0, column_comm);
}
/* if p is not a square number (grid_rows != grid_cols) */
else { } //not implemented
}
