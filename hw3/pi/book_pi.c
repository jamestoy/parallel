#include <stdio.h>
#include <time.h>
#include <assert.h>
#define INTERVALS 100000000

int
main() {
    clock_t start, stop;
    double t = 0.0;
    assert((start = clock()) != 1);

    double area, ysum, xi;
    int i;

    ysum = 0.0;
    for(i = 0; i < INTERVALS; i++) {
        xi = (1.0 / INTERVALS) * (i + 0.5);
        ysum += 4.0 / (1.0 + xi * xi);
    }
    area = ysum * (1.0 / INTERVALS);
    printf("area is %13.11f\n", area);

    stop = clock();
    t = (double)(stop-start)/CLOCKS_PER_SEC;
    printf("run time (in seconds): %f\n",t);

    return 0;
}
