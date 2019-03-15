#ifndef __L1C_MEMORY__
#define __L1C_MEMORY__

double* malloc_double(int N);
void free_double(double *x);

double** malloc_double_2D(l1c_int nrow, l1c_int ncol);

void free_double_2D(int nrow, double **dwork);

#endif
