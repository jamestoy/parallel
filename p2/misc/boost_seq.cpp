#include <iostream>

#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/timer.hpp>

namespace ublas = boost::numeric::ublas;

using namespace std;
//using namespace ublas;

void print_mat(ublas::matrix<double> mat, int row, int col)
{
    int i,j;
    for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			cout << " " << mat(i,j);
		}
		cout << endl;
	}
	cout << endl;
}

int main(int argc, char *argv[]) {

	int i,j,k;
	int l,m,n;

	if(argc != 2) {
                printf("Usage: mm_seq <mat_size>\n");
                return 1;
        } else
                l = m = n = atoi(argv[1]);

	ublas::matrix<double> A(n,n), B(n,n), C(n,n);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			A(i,j) = drand48();
			B(i,j) = drand48();
		}

	cout << "Size of matricies: " << n << "x" << n << endl;

	boost::timer timer;
	double time;

	//cout << "Printing original matrix 'A' BEFORE transpose" << endl;
	//print_mat(A,n,n);

	cout << "Transposing A before multiplication..." << endl;

	timer.restart();
	B = trans(B);
	time = timer.elapsed();

	cout << "Transpose time: " << time << " seconds" << endl;

	C = prod(A,B);
	time = timer.elapsed();

	cout << "Product computed.... total time: " << time << " seconds" << endl << endl;

	//print_mat(A,n,n);
	//print_mat(B,n,n);
	//print_mat(C,n,n);

	return 0;
}

