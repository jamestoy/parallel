#include <stdio.h>
#include <mpi.h>

#define TAG 3

int
main(int argc, char *argv[])
{
  int iter;
  int ar_size;
  if(argc != 3) {
    printf("Usage: commtest iterations array_size\n");
    return 1;
  } else {
    iter = atoi(argv[1]);
    ar_size = atoi(argv[2]);
  }

  int msg[ar_size]; int k;
  for(k = 0; k < ar_size; k++)
     msg[k] = 0;
  //printf("number of elements in array: %d\n",ar_size);
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
    MPI_Recv(&msg, ar_size, MPI_INT, from, TAG, MPI_COMM_WORLD, &stat);
    MPI_Send(&msg, ar_size, MPI_INT, to, TAG, MPI_COMM_WORLD);
  }
  } else {
  start = MPI_Wtime();
  for(j = 0; j < iter; j++) {
    MPI_Send(&msg, ar_size, MPI_INT, to, TAG, MPI_COMM_WORLD);
    MPI_Recv(&msg, ar_size, MPI_INT, from, TAG, MPI_COMM_WORLD, &stat);
  }
  end = MPI_Wtime();
  printf("run time: %f\n",end-start);
  }

  MPI_Finalize();
  return 0;
}
