#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 200

/* Prototypes */
double **mm(double **a, double **b, double **c, int l, int m, int n);
double **mm_recur(int crow, int ccol, 
				  int arow, int a col, 
				  int brow, int bcol,
				  int l,
				  int m,
				  int n);
double **init_mat(double **mat, int row, int col);
double **transpose(double **mat, int row, int col);
void print_mat(double **mat, int row, int col);

double **transpose(double **mat, int row, int col)
{
	int i,j;
	for (i = 0; i < row; i++) {
		for (j = i + 1; j < col; j++) {
			double tmp = mat[i][j];
			mat[i][j] = mat[j][i];
			mat[j][i] = tmp;
		}
	}
	return mat;
}

/* Computation step */
double **mm(double **a, double **b, double **c, int l, int m, int n)
{
	int i,j,k,sum;
	
//#pragma omp parallel for private(j,k)
	for (i = 0; i < l; i++) {
		for (j = 0; j < n; j++) {
			c[i][j] = 0;
//#pragma omp parallel for reduction(+:sum)
			for (k = 0; k < m; k++) {
				//sum += a[i][k] * b[k][j];
				c[i][j] += a[i][k] * b[k][j];
			}
		//c[i][j] = sum;
		}
	}
	return c;
}

double **mm_recur(int crow, int ccol, 
					 int arow, int a col, 
					 int brow, int bcol,
					 int l,
					 int m,
					 int n)
{
	int lhalf[3], mhalf[3], nhalf[3]; /* Quadrant sizes */
	int i, j, k; double *aptr, *bptr, *cptr;
	
	if (m * n > THRESHOLD) {
		
		/* B doesn't fit in cache---multiply block of A, B */
		
		lhalf[0] = 0; lhalf[1] = 1/2; lhalf[2] = l - l/2;
		mhalf[0] = 0; mhalf[1] = 1/2; mhalf[2] = m - m/2;
		nhalf[0] = 0; nhalf[1] = 1/2; nhalf[2] = n - n/2;
		for (i = 0; i < 2; i++)
			for (j = 0; j < 2; j++)
				for (k = 0; k < 2; k++)
					mm2 (crow+lhalf[i], ccol+mhalf[j],
					  arow+lhalf[i], acol+mhalf[k],
					  brow+mhalf[k], bcol+nhalf[j],
					  lhalf[i+1], mhalf[k+1], nhalf[j+1]);
	} else {
		
		/* B fits in cache---do standard multiply */	
		for (i = 0; i < 1; i ++)
//#pragma omp parallel for private(k)
			for (j = 0; j < n; j++) {
				cptr = &c[crow+i][ccol+j];
				aptr = &a[arow+i][acol];
				bptr = &b[brow][bcol+j];
				for (k = 0; k < m; k++) {
					*cptr += *(aptr++) * *bptr; bptr += N;
	return c;
}

void print_mat(double **mat, int row, int col)
{
    int i,j;
    for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			printf(" %f", mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

double **init_mat(double **mat, int row, int col)
{
	int i,j;
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			mat[i][j] = drand48();
		}
	}
	return mat;
}

int main(int argc, char *argv[])
{
	double **a, **b, **ret;
	clock_t start, stop;
	int i;
	int l, m, n;
	l = N;
	m = N;
	n = N;
	
	a = malloc(l * sizeof(double *));
	b = malloc(m * sizeof(double *));
	ret = malloc(l * sizeof(double *));
	if(!(a && b && ret)) {
        	printf("Out of memory; exiting...\n");
        	return 1;
    	}
	
	for (i = 0; i < l; i++)
		a[i] = malloc(m * sizeof(double));
	for (i = 0; i < m; i++)
		b[i] = malloc(n * sizeof(double));
	for (i = 0; i < l; i++)
		ret[i] = malloc(n * sizeof(double));
	
	a = init_mat(a,l,m);
	b = init_mat(b,m,n);
	print_mat(a,l,m);
	print_mat(b,m,n);
	
	if((start = clock()) == -1) {
		printf("Could not start the clock; exiting...\n");
		return 1;
	}
	
	ret = mm(a,b,ret,l,m,n);
	stop = clock();
	
	print_mat(ret,l,n);
	double time = (double)(stop - start)/CLOCKS_PER_SEC;
	printf("Elapsed time: %f\n", time);

	return 0;
}
