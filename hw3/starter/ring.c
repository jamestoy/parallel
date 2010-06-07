#include <stdio.h> 
#include <mpi.h>


int
main(int argc, char *argv[])
{
  int count = 5;
  int i = 0;
  int msg = 123;
  int rank, size;
  int from, to;
  MPI_Status stat;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  to   = (rank + 1) % size;
  from = (size + rank - 1) % size;
  
  printf("I am node %d of %d\n", rank, size);
  printf("Sending to %d and receiving from %d\n", to, from);
  
  if (rank == size - 1)
    MPI_Send(&msg, 1, MPI_INT, to, 4, MPI_COMM_WORLD);
  
  for (i = 0; i < count; i++) {
    MPI_Recv(&msg, 1, MPI_INT, from, 4, MPI_COMM_WORLD, &stat);
    printf("Node %d received %d\n", rank, msg);
    MPI_Send(&msg, 1, MPI_INT, to, 4, MPI_COMM_WORLD);
  }

  if (rank == 0) {
    MPI_Recv(&msg, 1, MPI_INT, from, 4, MPI_COMM_WORLD, &stat);
    printf("Node %d received %d\n", rank, msg);
  }
  MPI_Finalize();
  
  return 0;
}
