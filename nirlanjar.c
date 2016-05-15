#include "num.h"
#include <stdio.h>

#define PI 3.141592654

int main() {
    /**** SOAL 1 ****/
    printf("\nsoal 1\n");
    printf("a)\n");
    printf("bisection: %.5f\n", bisection(f1a, 0.001, 0.999));
    printf("secant: %.5f\n", secant(f1a, 0.001, 0.999));
    printf("false position: %.5f\n", falsepos(f1a, 0.001, 0.999));

    printf("\nb)\n");
    printf("bisection: %.5f\n", bisection(f1b, 0.001, 2.999*PI));
    printf("secant: %.5f\n", secant(f1b, 0.001, 2.999*PI));
    printf("false position: %.5f\n", falsepos(f1b, 0.001, 2.999*PI));

    printf("\nc)\n");
    printf("bisection: %.5f\n", bisection(f1c, 0.001, 0.999*PI));
    printf("secant: %.5f\n", secant(f1c, 0.001, 0.999*PI));
    printf("false position: %.5f\n", falsepos(f1c, 0.001, 0.999*PI));

    /**** SOAL 2 ****/
    printf("\nsoal 2\n");
    printf("newton raphson: %.5f\n", newton(f2, df2, 0.01));
    printf("bisection: %.5f\n", bisection(f2, 0.001, 100));
    printf("secant: %.5f\n", secant(f2, 0.001, 100));
    printf("false position: %.5f\n", falsepos(f2, 0.001, 100));

    /**** SOAL 3 ****/
    printf("\nsoal 3\n");
    printf("newton raphson: %.5f\n", newton(f3, df3, 0.01));
    printf("bisection: %.5f\n", bisection(f3, 0.001, 100));
    printf("secant: %.5f\n", secant(f3, 0.001, 100));
    printf("false position: %.5f\n", falsepos(f3, 0.001, 100));

    return 0;
}
