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
/* 
// Based on loop nest optimization
// we could also add blocking to increase cache performance
// which would ensure we are only putting slices in that are the size
// of the cache......
double **mm(double **a, double **b, double **c, int l, int m, int n)
{

	int acc00, acc01, acc10, acc11;
	int i,j,k
	for (i = 0; i < l; i += 2)
	{
	    for (j = 0; j < n; j += 2)
	    {
	        acc00 = acc01 = acc10 = acc11 = 0;
	        for (k = 0; k < m; k++)
	        {
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
/* Computation step */
double **mm(double **a, double **b, double **c, int l, int m, int n)
{
	int i,j,k;
	
	for (i = 0; i < l; i+=3) {
		for (j = 0; j < n; j+=3) {
			c[i][j] = 0;
			for (k = 0; k < m; k+=3) {
				c[i][j] += a[i][k] * b[k][j];
			}
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
	} else
		l = m = n = atoi(argv[1]);
	
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
	
	printf("matrix size: %d x %d\n", n,n);
	
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
