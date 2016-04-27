#ifndef __NUM_H__
#define __NUM_H__

#include "matrix.h"

void pivot(matrix *m, long cur_row);
void gaussjordan(matrix *a, double *b, double *results);
void gaussnaive(matrix *a, double *b, double *results);
void doolittle(matrix *a, double *b, double *results);
void crout(matrix *a, double *b, double *results);

typedef double (*function)(double x);
double f1a(double);
double f1b(double);
double f1c(double);
double f2(double);
double df2(double);
double f3(double);
double df3(double);
double bisection(function, double, double);
double falsepos(function, double, double);
double secant(function, double, double);
double newton(function, function, double);

#endif
