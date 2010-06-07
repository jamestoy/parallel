#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <errno.h>

int
main(int argc, char *argv[])
{
  int msg;
  if(argc == 1) {
    printf("please input exactly one arguement (an integer) to set msg size\n");
    exit(1);
  } else {
    msg = atoi(argv[1]);
  }
  int count = 5;
  int i = 0;
  int rank, size;
  int from, to;
  MPI_Status stat;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  to   = (rank + 1) % size;
  from = (size + rank - 1) % size;

  char hostname[128];
  if(gethostname(hostname, sizeof hostname) < 0)
    printf("gethostname errno: %s\n",strerror(errno));

  printf("node name: %s\n", hostname);
  printf("node %d of %d\n", rank, size);
  printf("Sending to %d and receiving from %d\n", to, from);
  
  if (rank == size - 1) {
    msg++;	//init incrememnt
    MPI_Send(&msg, 1, MPI_INT, to, 4, MPI_COMM_WORLD);
  }

  for (i = 0; i < count; i++) {
    MPI_Recv(&msg, 1, MPI_INT, from, 4, MPI_COMM_WORLD, &stat);
    printf("Node %d received %d\n", rank, msg);
    msg++;	//increment (part c)
    MPI_Send(&msg, 1, MPI_INT, to, 4, MPI_COMM_WORLD);
  }

  if (rank == 0) {
    MPI_Recv(&msg, 1, MPI_INT, from, 4, MPI_COMM_WORLD, &stat);
    printf("Node %d received %d\n", rank, msg);
  }
  printf("node %s is all finished\n", hostname);
  MPI_Finalize();
  
  return 0;
}
