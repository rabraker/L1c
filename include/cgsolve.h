#ifndef _CGSOLVE_
#define _CGSOLVE_


typedef struct CgResults_{
  double cgres;
  int cgiter;
} CgResults;

typedef struct CgParams_{
  int verbose;
  int max_iter;
  double cgres;
  double tol;
} CgParams;


int cgsolve(double *x, double *b, size_t n_b,  double *Dwork,
            void(*AX_func)(int n, double *x, double *b, void *AX_data),
            void *AX_data, CgResults *cg_result, CgParams cg_params);

void dgemv_RowOrder(double *A, int m_A, int n_A, double *x, double *b);

void Ax(int n, double *x, double *b, void *AX_data);



#endif
