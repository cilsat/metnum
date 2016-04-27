#include "num.h"
#include <stdio.h>

#define PI 3.141592654

int main() {
    /**** SOAL 1 ****/
    printf("\nsoal 1\n");
    printf("bisection: %.5f\n", bisection(f1a, 0.001, 1));
    printf("secant: %.5f\n", secant(f1a, 0.001, 1));
    printf("false position: %.5f\n", falsepos(f1a, 0.001, 1));

    printf("bisection: %.5f\n", bisection(f1b, 0.001, 3*PI));
    printf("secant: %.5f\n", secant(f1b, 0.001, 3*PI));
    printf("false position: %.5f\n", falsepos(f1b, 0.001, 3*PI));

    printf("bisection: %.5f\n", bisection(f1c, 0.001, 10));
    printf("secant: %.5f\n", secant(f1c, 0.001, 10));
    printf("false position: %.5f\n", falsepos(f1c, 0.001, 10));

    /**** SOAL 2 ****/
    printf("\nsoal 2\n");
    printf("newton raphson: %.5f\n", newton(f2, df2, 0.01));

    /**** SOAL 3 ****/
    printf("\nsoal 3\n");
    printf("newton raphson: %.5f\n", newton(f3, df3, 0.01));

    return 0;
}
