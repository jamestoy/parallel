/****************************************/
/* Compute pi to arbitrary precision    */
/* Author Roy Williams February 1994    */
/* Uses Machin's formula...             */
/* pi/4 = 4 arctan(1/5) - arctan(1/239) */
/****************************************/
/* compile with cc -O -o pi pi.c        */
/* run as "pi 1000"                     */
/****************************************/
/* The last few digits may be wrong.......... */

#include <stdio.h>

#define BASE    10000
int nblock;
int *tot;
int *t1;
int *t2;
int *t3;

main(argc, argv)
int argc;
char **argv;
{
        int ndigit = 60;
        if(argc == 2)
                ndigit = atoi(argv[1]);
        else {
                fprintf(stderr, "Usage: %s ndigit\n", argv[0]);
                fprintf(stderr, "(Assuming 60 digits)\n");
        }
        if(ndigit < 20) ndigit = 20;
        nblock = ndigit/4;
        tot = (int *)malloc(nblock*sizeof(int));
        t1 = (int *)malloc(nblock*sizeof(int));
        t2 = (int *)malloc(nblock*sizeof(int));
        t3 = (int *)malloc(nblock*sizeof(int));
        if(!tot || !t1 || !t2 || !t3){
                fprintf(stderr, "Not enough memory\n");
                exit(1);
        }

        arctan(tot, t1, t2, 5, 1);
        mult(tot, 4);
        arctan(t3, t1, t2, 239, 2);
        sub(tot, t3);
        mult(tot, 4);
        print(tot);
}

arctan(result, w1, w2, denom, onestep)
int *result, *w1, *w2, denom, onestep;
{
        int denom2 = denom*denom;
        int k = 1;

        set(result, 1);
        div(result, denom);
        copy(w1, result);

        do{
                if(onestep)
                        div(w1, denom2);
                else {
                        div(w1, denom);
                        div(w1, denom);
                }
                copy(w2, w1);
                div(w2, 2*k+1);
                if(k%2)
                        sub(result, w2);
                else
                        add(result, w2);
                k++;
        } while(!zero(w2));
}

copy(result, from)
int *result, *from;
{
        int i;
        for(i=0; i<nblock; i++)
                result[i] = from[i];
}

zero(result)
int *result;
{
        int i;
        for(i=0; i<nblock; i++)
                if(result[i])
                        return 0;
        return 1;
}

add(result, increm)
int *result, *increm;
{
        int i;
        for(i=nblock-1; i>=0; i--){
                result[i] += increm[i];
                if(result[i] >= BASE){
                        result[i] -= BASE;
                        result[i-1]++;
                }
        }
}

sub(result, decrem)
int *result, *decrem;
{
        int i;
        for(i=nblock-1; i>=0; i--){
                result[i] -= decrem[i];

                if(result[i] < 0){
                        result[i] += BASE;
                        result[i-1]--;
                }
        }
}

mult(result, factor)
int *result, factor;
{
        int i, carry = 0;
        for(i=nblock-1; i>=0; i--){
                result[i] *= factor;
                result[i] += carry;
                carry = result[i]/BASE;
                result[i] %= BASE;
        }
}

div(result, denom)
int *result, denom;
{
        int i, carry = 0;

        for(i=0; i<nblock; i++){
                result[i] += carry*BASE;
                carry = result[i] % denom;
                result[i] /= denom;
        }
}

set(result, rhs)
int *result, rhs;
{
        int i;
        for(i=0; i<nblock; i++)
                result[i] = 0;
        result[0] = rhs;
}

print(result)
int *result;
{
        int i, k;
        char s[10];
        printf("%1d.\n", result[0]);
        for(i=1; i<nblock; i++){
                sprintf(s, "%4d ", result[i]);
                for(k=0; k<5; k++)
                        if(s[k] == ' ') s[k] = '0';
                printf("%c%c%c%c", s[0], s[1], s[2], s[3]);
                if(i%15 == 0) printf("\n");
        }
        printf("\n");
}

