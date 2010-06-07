//#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static int n;	 	 /* matrix size */
static int bsize;	 /* submatrix size */
static int bnum;	 /* per row/col submatricies in main matrix */

static int nthread;	 /* number of threads */
static int sthread;      /* set threads to this number */
static int myid; 	 /* whoami */

static double	*A,	 /* matricies */
		*B,
		*C,
		*tmp; 

/* map for non block version */
int LOCAL(int i)
{
	return i%(n/nthread);
}

/* map for block version */
int LOCALB(int i)
{
	return i%(bnum/nthread);
}

/* get pointer to start of matricies */
double *MAP_BLK_A(int i, int j)
{	
	return A+(i*bnum+j)*bsize*bsize;
}

double *MAP_BLK_B(int i, int j)
{	
	return B+(j*bnum+i)*bsize*bsize;
}

double *MAP_BLK_C(int i, int j)
{	
	return C+(i*bnum+j)*bsize*bsize;
}

double *MAP_BLK_BUF(int k)
{	
	return tmp+k*bsize*bsize;
}

/* owners of a block */
int OWNER_A(int i, int j)
{
	return i/(n/nthread);
}

int OWNER_B(int i, int j)
{
	return j/(n/nthread);
}

int OWNER_C(int i, int j)
{
	return i/(n/nthread);
}

/* get the addresses of A/B/C */
double * MAP_ELEM_A(int k, int l)
{
	int local_i=LOCAL(k), local_j=l;
	int local_bi=local_i/bsize, local_bj=local_j/bsize;
	int blk_i=local_i%bsize, blk_j=local_j%bsize;
	return MAP_BLK_A(local_bi, local_bj)+blk_i*bsize+blk_j;
}

double * MAP_ELEM_B(int k, int l)
{
	int local_i=k, local_j=LOCAL(l);
	int local_bi=local_i/bsize, local_bj=local_j/bsize;
	int blk_i=local_i%bsize, blk_j=local_j%bsize;
	return MAP_BLK_B(local_bi, local_bj)+blk_i*bsize+blk_j;
}

double * MAP_ELEM_C(int k, int l)
{
	int local_i=LOCAL(k), local_j=l;
	int local_bi=local_i/bsize, local_bj=local_j/bsize;
	int blk_i=local_i%bsize, blk_j=local_j%bsize;
	return MAP_BLK_C(local_bi, local_bj)+blk_i*bsize+blk_j;
}

/* Based on loop nest optimization -- we could also add blocking to increase cache
performance
double **mm(double *c, double *a, double *b)
{

int acc00, acc01, acc10, acc11;
int i,j,k
for (i = 0; i < bsize; i += 2)
{
    for (j = 0; j < bsize; j += 2)
    {
        acc00 = acc01 = acc10 = acc11 = 0;
        for (k = 0; k < bsize; k++)
        {
	   //c[i*bsize+j]+=a[i*bsize+k]*b[k*bsize+j];
            acc00 += a[k][j + 0] * b[i + 0][k];
            acc01 += a[k][j + 1] * b[i + 0][k];
            acc10 += a[k][j + 0] * b[i + 1][k];
            acc11 += a[k][j + 1] * b[i + 1][k];
        }
        c[i + 0][j + 0] = acc00;
        c[i + 0][j + 1] = acc01;
        c[i + 1][j + 0] = acc10;
        c[i + 1][j + 1] = acc11;
    }
}
}
*/

void do_sub_mm(double *c, double *a, double *b)
{
	int i, j, k;

#pragma omp parallel for private(j,k)
	for (i=0; i<bsize; i++)
		for (j=0; j<bsize; j++)
			for (k=0; k<bsize; k++)
				c[i*bsize+j]+=a[i*bsize+k]*b[k*bsize+j];
}

void print_result()
{
	int i, j;

	for (i=0; i<n; i++)
		if (OWNER_C(i, 0)==myid) {
			for (j=0; j<n; j++)
				printf("%4f  ", *MAP_ELEM_C(i,j));
			printf("\n");
		}
}

int main(int argc, char *argv[])
{
	int i, j, k, token;
	clock_t start, stop;
//	//MPI_Status status;
//	//MPI_Init(&argc, &argv);
//	//MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//	//MPI_Comm_size(MPI_COMM_WORLD, &nthread);
	myid = omp_get_thread_num();

	if (argc != 4) {
		if (myid==0)
			fprintf(stderr, "Usage, %s <n> <n/q> <num_threads>.\n", argv[0]);
//		//MPI_Finalize();
		return -1;
	}

	n=atoi(argv[1]);
	bsize=atoi(argv[2]);
	sthread = atoi(argv[3]);

	omp_set_num_threads(sthread);	
	nthread = omp_get_num_threads();
	
	printf("matrix size: %dx%d ... mesh size: %dx%d ... threads: %d\n",n,n,bsize,bsize,sthread);
	
	bnum=n/bsize;
	if ((bnum*bsize!=n) || (bnum/nthread*nthread!=bnum)) {
		if (myid==0)
			fprintf(stderr, "Invalid parameters <nthread=%d> <n=%d> <n/q=%d>.\n", nthread, n, bsize);
//		//MPI_Finalize();
		return -1;
	}

	A=(double *)malloc(sizeof(double)*(n*n/nthread));
	B=(double *)malloc(sizeof(double)*(n*n/nthread));
	C=(double *)malloc(sizeof(double)*(n*n/nthread));
	tmp=(double *)malloc(sizeof(double)*(bsize*n));

	if (!(A && B && C && tmp)) {
		if (myid==0)
			fprintf(stderr, "Out of Memory.\n");
		//MPI_Finalize();
		return -1;
	}

	memset((char *)C, 0, sizeof(double)*(n*n/nthread));

	if (myid==0)
		printf("Initialization ...\n");
	

	/* initialize A and B*/
	for (i=0; i<n; i++)
		for (j=0; j<n; j++) {
			if (OWNER_A(i,j)==myid) {
				*MAP_ELEM_A(i,j)=j;
			}
			if (OWNER_B(i,j)==myid) {
				*MAP_ELEM_B(i,j)=i+j;
			}
	}

	if (myid==0) {
		printf("Done!\n\n");
		printf("Multiplication ...\n");
	}

	if ((start = clock()) == -1) {
		printf("Unable to start the clock() -- exiting...\n");
		return 1;
	}

	//double mpi_start = MPI_Wtime();

#pragma omp parallel for private(i,k)
	for (j=0; j<bnum; j++)
	{
		myid = omp_get_thread_num();
		//printf("after pragma myid = %d\n",myid);
		int root=OWNER_B(0, j*bsize);
		if (root==myid)
			memcpy(tmp, MAP_BLK_B(0, LOCALB(j)), sizeof(double)*(bsize*n));

		//MPI_Bcast(tmp, bsize*n, MPI_DOUBLE, root, MPI_COMM_WORLD);

		for (i=0; i<bnum/nthread; i++) {
			/* update C[i][j] */
			for (k=0; k<bnum; k++) 
				/* C[i][j]+=A[i][k]*tmp[k] */
				do_sub_mm(MAP_BLK_C(i, j), MAP_BLK_A(i, k), MAP_BLK_BUF(k));
		}
	}
	stop = clock();
	//double mpi_stop = MPI_Wtime();
	//double mpi_time = stop = start;
	//printf("MPI time : %f\n",mpi_time);
	double time = ((double)(stop-start)/CLOCKS_PER_SEC)/sthread;
	printf("time = %f\n",time);
	if (myid==0)
		printf("Done!\n\n");

	/* If the size of the matrix is greater than 10 */
//	if (n<=10) {
//		if (myid==0)
//			printf("Printing Matrix C ...\n");

		/* implementing a token pass protocol to print the C */
//		if (myid!=0) {
//			MPI_Recv(&token, 1, MPI_INT, myid-1, 2, MPI_COMM_WORLD, &status);
//		
	
	//	print_result();
		fflush(stdout);

//		if (myid<nthread-1)
//			MPI_Send(&token, 1, MPI_INT, myid+1, 2, MPI_COMM_WORLD);
//		if (myid==nthread-1)
//			printf("Done!\n");
//	}

	//MPI_Finalize();

	return 0;
}
