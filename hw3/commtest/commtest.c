#include <stdio.h>
#include <mpi.h>

#define NUM_ELMNTS	1
#define TAG		3

int
main(int argc, char *argv[])
{
  int iter = 0;
  if(argc != 2) {
    printf("Usage: commtest <integer>\n");
    return 1;
  } else {
    iter = atoi(argv[1]);
  }

  int msg = 33;
  int i, j, rank, size, from, to;
  double start, end;
  MPI_Status stat;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  to   = (rank + 1) % size;
  from = (size + rank - 1) % size;

  printf("I am node %d of %d\n", rank+1, size);
  printf("Sending to %d and receiving from %d\n", to, from);  

  if(rank == 0) {
  for (i = 0; i < iter; i++) {
    MPI_Recv(&msg, 1, MPI_INT, from, TAG, MPI_COMM_WORLD, &stat);
    //printf("Node %d received %d\n", rank, msg);
    MPI_Send(&msg, 1, MPI_INT, to, TAG, MPI_COMM_WORLD);
  }
  } else {
  start = MPI_Wtime();
  for(j = 0; j < iter; j++) {
    MPI_Send(&msg, NUM_ELMNTS, MPI_INT, to, TAG, MPI_COMM_WORLD);
    MPI_Recv(&msg, NUM_ELMNTS, MPI_INT, from, TAG, MPI_COMM_WORLD, &stat);
  }
  end = MPI_Wtime();
  printf("run time: %f\n",end-start);
  }

  MPI_Finalize();
  return 0;
}
