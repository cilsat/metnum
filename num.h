#ifndef __NUM_H__
#define __NUM_H__

#include "matrix.h"

void pivot(matrix *m, long cur_row);
void gaussjordan(matrix *a, double *b, double *results);
void gaussnaive(matrix *a, double *b, double *results);
void ludecomp(matrix *a, double *b, double *results);

typedef double (*nirlanjar)(double x);
double f1a(double);
double f1b(double);
double f1c(double);
double f2(double);
double df2(double);
double f3(double);
double df3(double);
double bisection(nirlanjar, double, double);
double falsepos(nirlanjar, double, double);
double secant(nirlanjar, double, double);
double newton(nirlanjar, nirlanjar, double);

#endif
