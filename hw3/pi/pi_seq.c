#include <stdio.h>
#include <time.h>
#include <assert.h>

int main() {
    clock_t start, stop;
    double t = 0.0;
    assert((start = clock()) != -1);

    int r[2800 + 1];
    int i, k;
    int b, d;
    int c = 0;

    for (i = 0; i < 2800; i++) {
        r[i] = 2000;
    }

    for (k = 2800; k > 0; k -= 14) {
        d = 0;

        i = k;
        for (;;) {
            d += r[i] * 10000;
            b = 2 * i - 1;

            r[i] = d % b;
            d /= b;
            i--;
            if (i == 0) break;
            d *= i;
        }
        printf("%.4d", c + d / 10000);
        c = d % 10000;
    }

    stop = clock();
    t = (double)(stop-start)/CLOCKS_PER_SEC;
    printf("\n\nrun time (in seconds): %f\n",t);

    return 0;
}
