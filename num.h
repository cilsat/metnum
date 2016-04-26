#ifndef __NUM_H__
#define __NUM_H__

#include "matrix.h"

void gaussjordan(matrix *a, float *b, float *results);
void gaussnaive(matrix *a, float *b, float *results);
void ludecomp(matrix *a, float *b, float *results);

#endif
