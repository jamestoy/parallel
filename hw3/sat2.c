/*
 *   Circuit Satisfiability, Version 2
 *
 *   This enhanced version of the program prints the
 *   total number of solutions.
 */

#include <mpi.h>
#include <stdio.h>
#define NUM_VARS 23
int
main (int argc, char *argv[]) {
   int count;            /* Solutions found by this proc */
   int global_count;     /* Total number of solutions */
   int i;
   int id;               /* Process rank */
   int p;                /* Number of processes */
   int check_circuit (int, int);
   double start, end, time, global_time;

   MPI_Init (&argc, &argv);
   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);

   count = 0;
   time = 0.0;
   global_time = 0.0;

   start = MPI_Wtime();
   for (i = id; i < (1 << NUM_VARS); i += p)
      count += check_circuit (id, i);
   end = MPI_Wtime();
   time += (end - start);

   printf("process %d took %f time to do the main loop\n", id, time);
   MPI_Reduce(&time, &global_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

//COUNTS
   MPI_Barrier(MPI_COMM_WORLD); /* for clarity's sake */
   printf("process %d has %d solutions\n", id, count);
   //MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

   MPI_Barrier(MPI_COMM_WORLD);
   printf ("process %d is done\n", id);
   fflush (stdout);
   MPI_Finalize();
   if(!id) printf("total time: %f\n", global_time);
   /*if (!id) printf ("There are %d different solutions\n",global_count);*/

   return 0;
}

/* Return 1 if 'i'th bit of 'n' is 1; 0 otherwise */
#define EXTRACT_BIT(n,i) ((n&(1<<i))?1:0)

int check_circuit (int id, int z) {
   int v[NUM_VARS];        /* Each element is a bit of z */
   int i;

   for (i = 0; i < NUM_VARS; i++) v[i] = EXTRACT_BIT(z,i);
   if ((v[0] || v[1]) && (!v[1] || !v[3]) && (v[2] || v[3])
      && (!v[3] || !v[4]) && (v[4] || !v[5])
      && (v[5] || !v[6]) && (v[5] || v[6])
      && (v[6] || !v[15]) && (v[7] || !v[8])
      && (!v[7] || !v[13]) && (v[8] || v[9])
      && (v[8] || !v[9]) && (!v[9] || !v[10])
      && (v[9] || v[11]) && (v[10] || v[11])
      && (v[12] || v[13]) && (v[13] || !v[14])
      && (v[14] || v[15]) && (v[15] || !v[16])
      && (v[16] || v[17]) && (v[17] || !v[18])
      && (v[17] || v[20]) && (v[18] || !v[19])
      && (v[19] || v[20]) && (v[21] || !v[22])
      && (v[22] || v[23]) && (v[20] || !v[21])) {
      printf ("%d) %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d\n", id,
         v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],
         v[10],v[11],v[12],v[13],v[14],v[15],v[16],v[17],v[18],v[19],v[20],v[21],v[22],v[23]);
      fflush (stdout);
      return 1;
   } else return 0;
}
