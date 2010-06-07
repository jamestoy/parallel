#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <omp.h>

#define N 32 

/* PROTOS */
void a(double aa[], double ba[], double x);
void b(double aa[], double ba[], int n);
void c(double aa[], int n);
void d(double aa[], double ba[], int n);
void e(double aa[], double ba[], int n);
void f(double aa[], double ba[], int n);
void g(double aa[], int k);
void h(double aa[], double ba[], double b, int k, int n);
double foo(int i);

void a(double aa[], double ba[], double x)
{
	int i,n;
    /* checks to make sure we are in bounds
     * this is just like all the others seen below
     */
	if(((int) sqrt(x)) < N) n = (int) sqrt(x);
	else n = N;

    /* This is a ||izable for loop for sure.
     * all it is really doing is doing an assignment.
     * Additionally the sqrt was broken out and a
     * test was added to ensure no overflow.  Other
     * than that this seems a very simple example.
     * OpenMP will schedule iterations to threads, 
     * and presto! ||ized.
     */
	#pragma omp parallel for
	for(i = 0; i < n; i++) {
		aa[i] = 2.3 * i;
		if(i < 10) ba[i] = aa[i];
	}
}

void b(double aa[], double ba[], int n)
{
	int flag = 0;
	int i;

	if(n > 32) n = N;

	/* This is not parallelizable because of the fact
	 * that a flag is set based on a conditional statement.
	 * The flag is presumably to mark an important action; and,
	 * if parallelized, there is the potential for several threads
	 * to be in the same position and this can lead to unpredictable
	 * behavior depending upon what the flag is for (termination, 
	 * start a function, etc). So this is why I consider this
	 * for loop unparallelizable.
 	 */ 
	for (i = 0; (i < n) & (!flag); i++) {
		aa[i] = 2.3 * i;
		if( aa[i] < ba[i]) flag = 1;
	}
}

void c(double aa[], int n)
{
	int i;
	if(n > N) n = N;
    /* This is parallelizable because it is a simple
     * assignment statement from a function call.
     * Since A) each iteration is given to a thread
     * and   B) each call to foo is based on this
     * we know that we can parallelize this safely
     * with the parallel for compiler directive.
     * Also, though a bit unfounded perhaps, we put
     * in schedule(guided) since we are not sure
     * how long foo could take or what it does
     * this will dynamically allocate with self-checking
     * heuristics to ensure speed with a min chunk size of 1
     */
	#pragma omp parallel for schedule(guided)
	for (i = 0; i < n; i++)
		aa[i] = foo(i);
}

void d(double aa[], double ba[], int n)
{
	int i;
	if(n > N) n = N;
    /* This is ||izable because again this is an
     * assignment statement; the only weirdness
     * about this function is the fact we do not know what
     * the function "foo" is doing and so we utilize the 
     * schedule(guided) to ensure if foo is a crazy computation
     * we are being careful about how threads are allocated.
     */
	#pragma omp parallel for schedule(guided)
	for (i = 0; i < n; i++)
	{
		aa[i] = foo(i);
		if(aa[i] < ba[i]) aa[i] = ba[i];
	}
}

void e(double aa[], double ba[], int n)
{
	int i;
	if(n > N) n = N;

	/* This is a more specific version of function b.
	 * The reason for this is because of the fact that
	 * the termination of the function is contingent upon
	 * a conditional that cannot be parallelized due 
	 * to the fact that termination is conditional. 
	 * The addition of parallelism to this function could
	 * easily lead to unpredictable / delayed behavior (for example
	 * if the break function was to be the only one delayed in
	 * a context switch -- i.e. the condition to stop had already
	 * been met but that thread was not running because it was 
	 * waiting....)
	 */
	for(i = 0; i < n; i++)
	{
		aa[i] = foo(i);
		if (aa[i] < ba[i]) break;
	}
}

void f(double aa[], double ba[], int n)
{
    /* This function is parallelizable AND
     * we get to make use of the reduction clause
     * in OpenMP.  We could put a critical section
     * around the dotp but we can see that using
     * a reduction clause is far superior in all
     * cases since it is not only faster on a single
     * thread but when additional threads are added we
     * get speedup. Also the reduction is necessary because
     * without it we can get nondeterministic behavior, namely
     * two threads read and compute the same sum before the first
     * one has written back the value.  This will produce failure.
     * Finally note the if clause in the pragma, this will notify
     * the compiler whether or not the loop should be
     * executed in parallel or sequentially.  This is useful
     * because there are often times when the time taken
     * to do the forking and joining is much greater than
     * the benefit of parallelizing.
     */
	int i;
	if(n > N) n = N;
	int dotp = 0;
	#pragma omp parallel for reduction(+:dotp) if(n > 5000)
	for (i = 0; i < n; i++)
		dotp += aa[i] * ba[i];
}

void g(double aa[], int k)
{
    /* This is parallelizable because of the fact
     * that all that is being done is an assignment.
     * The threads will be assigned to the simple assignment
     * of an array indecie plus the array indicie minus a constant
     * there is no reason this cannot be done in parallel.
     */
	int i;
	if((2 * k) > N) k = N/2;
	#pragma omp parallel for
	for(i = k; i < 2*k; i++)
		aa[i] = aa[i] + aa[i-k];
}

void h(double aa[], double ba[], double b, int k, int n)
{
	/* Not parallelizable because there is a possibility
     * that there is a data dependency problem.
     * The reference to aa[i-k] can be contingent
     * upon iterations that may not have been completed
     * as a result this does not work for parallization.
     */
    int i;
	if(n > N) n = N;
	for(i = k; i < n; i ++)
		aa[i] = b * aa[i-k];
}

/* Simple foo function to take a simple computation into account */
double foo(int i)
{
	return (i * 200.3 + 3.2 - 2.2);
}

/* Main with timings and array initialization */
int main(int argc, char *argv[])
{
	clock_t sa,ea,sb,eb,sc,ec,sd,ed,se,ee,sf,ef,sg,eg,sh,eh;
	double ta,tb,tc,td,te,tf,tg,th;
	int glob_i;

	//XXX: fix this.
	//unsigned short *oldptr;
	//unsigned short _seed[3];
	//unsigned short *seed48(unsigned short *);
	//oldptr = seed48(_seed);
	
	double aa[N]; double ba[N];
	for(glob_i = 0; glob_i < N; glob_i++) {
		//XXX: fix this.
		//aa[glob_i] = drand48(); ba[glob_i] = drand48();
		aa[glob_i] = ba[glob_i] = glob_i;
	}

	/*for(glob_i = 0; glob_i < N; glob_i++)
		printf("aa[i] = %f AND ba[i] = %f\n", aa[glob_i], ba[glob_i]);*/

	printf("OMP_NUM_THREADS = %d\n", omp_get_num_threads());
	
	assert((sa = clock()) != -1);
	a(aa,ba,N);
	ea = clock(); 
	ta = (double)(ea-sa)/CLOCKS_PER_SEC;
	printf("function a returns time: %f\n",ta);

	assert((sb = clock()) != -1);
        b(aa,ba,N);
        eb = clock();
	tb = (double)(eb-sb)/CLOCKS_PER_SEC;
	printf("function b returns time: %f\n",tb);

	assert((sc = clock()) != -1);
        c(aa,N);
        ec = clock();
	tc = (double)(ec-sc)/CLOCKS_PER_SEC;
	printf("function c returns time: %f\n",tc);

	assert((sd = clock()) != -1);
        d(aa,ba,N);
        ed = clock();
	td = (double)(ed-sd)/CLOCKS_PER_SEC;
	printf("function d returns time: %f\n",td);

	assert((se = clock()) != -1);
        e(aa,ba,N);
        ee = clock();
	te = (double)(ee-se)/CLOCKS_PER_SEC;
	printf("function e returns time: %f\n",te);

	assert((sf = clock()) != -1);
        f(aa,ba,N);
        ef = clock();
	tf = (double)(ef-sf)/CLOCKS_PER_SEC;
	printf("function f returns time: %f\n",tf);

	assert((sg = clock()) != -1);
        g(aa,N);
        eg = clock();
	tg = (double)(eg-sg)/CLOCKS_PER_SEC;
	printf("function g returns time: %f\n",tg);

	assert((sh = clock()) != -1);
        h(aa,ba,3.312312,2,N);
        eh = clock();
	th = (double)(eh-sh)/CLOCKS_PER_SEC;
	printf("function h returns time: %f\n",th);

	return 0;
}
