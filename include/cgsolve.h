#ifndef _CGSOLVE_
#define _CGSOLVE_

int cgsolve(double *x, double *b, size_t n_b, double tol,
                   int max_iter, int verbose, double *Dwork,
                   void(*AX_func)(int n, double *x, double *b, void *AX_data), void *AX_data);

void dgemv_RowOrder(double *A, int m_A, int n_A, double *x, double *b);

void Ax(int n, double *x, double *b, void *AX_data);

#endif
