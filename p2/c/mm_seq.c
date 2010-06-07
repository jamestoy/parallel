#include <stdio.h>
#include <stdlib.h>
#include <time.h>


/* Prototypes */
double **mm(double **a, double **b, double **c, int l, int m, int n);
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
	int i,j,k;
	double sum;	
//#pragma omp parallel for private(j,k)
	for (i = 0; i < l; i++) {
		for (j = 0; j < n; j++) {
			sum = 0;
//#pragma omp parallel for reduction(+:sum)
			for (k = 0; k < m; k++) {
				sum += a[i][k] * b[k][j];
			}
			c[i][j] = sum;
		}
	}
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

	if(argc != 2) {
                printf("Usage: mm_seq <mat_size>\n");
                return 1;
        } else {
                l = m = n = atoi(argv[1]);
                printf("Preparing matrix multiplication of an %dx%d matrix\n",n,n);
        }
	
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
	//print_mat(a,l,m);
	//print_mat(b,m,n);
	
	if((start = clock()) == -1) {
		printf("Could not start the clock; exiting...\n");
		return 1;
	}
	
	b = transpose(b,m,n);
	ret = mm(a,b,ret,l,m,n);
	stop = clock();
	
	//print_mat(ret,l,n);
	double time = (double)(stop - start)/CLOCKS_PER_SEC;
	printf("Elapsed time: %f\n", time);

	return 0;
}
