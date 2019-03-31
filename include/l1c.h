#ifndef __L1C__
#define __L1C__
// #include "config.h"
typedef int l1c_int;

#define DALIGN  64

double* l1c_malloc_double(int N);

void l1c_free_double(double *);

double** l1c_malloc_double_2D(l1c_int nrow, l1c_int ncol);

void l1c_free_double_2D(int nrow, double **dwork);


typedef struct l1c_AxFuns {
  /*Everybody must implement these.*/
  void(*Ax)(double *x, double *y);
  void(*Aty)(double *y, double *x);
  void(*AtAx)(double *x, double *z);

  /* Optional. Used by dct and dct2.*/
  void(*M)(double *x);
  void(*Mt)(double *y);

  /* Currently unused.*/
  void(*E)(double *x);
  void(*Et)(double *y);

  /* All implementations must define destroy*/
  void(*destroy)(void);

  /* Reserved for use.*/
  void *data;
}l1c_AxFuns;


/** Struct containing artifacts of the cgsolve routine. */
typedef struct l1c_CgResults_{
  /** Residual */
  double cgres;
  /** Number of completed conjugate gradient iterations. */
  l1c_int cgiter;

} l1c_CgResults;


typedef struct l1c_CgParams_{
  /** If 0, print nothing, if >0, print status every verbose-th iteration. */
  l1c_int verbose;
  /** Maximum number of solver iterations.*/
  l1c_int max_iter;
  /** Solver tolerance.*/
  double tol;
} l1c_CgParams;


typedef struct l1c_LBResult{
  double l1;
  int    total_newton_iter;
  int    total_cg_iter;
  int    status;

}l1c_LBResult;


typedef struct l1c_L1qcOpts{
  double epsilon;
  double mu;
  double lbtol;
  double tau;
  int lbiter;
  double newton_tol;
  int newton_max_iter;
  int verbose;
  double l1_tol;
  double cg_tol;
  int cg_maxiter;
  int cg_verbose;
  int warm_start_cg;
}l1c_L1qcOpts;

int l1c_cgsolve(l1c_int N, double *x, double *b, double **Dwork,
                void(*AX_func)(l1c_int n, double *x, double *b, void *AX_data),
                void *AX_data, l1c_CgResults *cg_result, l1c_CgParams cg_params);

int l1c_cgsolve_diag_precond(l1c_int N, double *x, double *b, double *M_inv_diag, double **Dwork,
                             void(*AX_func)(l1c_int n, double *x, double *b, void *AX_data), void *AX_data,
                             l1c_CgResults *cg_result, l1c_CgParams cg_params);

l1c_LBResult l1c_l1qc_newton(l1c_int N, double *x, l1c_int M, double *b,
                             l1c_L1qcOpts params, l1c_AxFuns Ax_funs);

int l1qc_dct(int Nrow, int Ncol, double *x_out, int M, double *b, l1c_int *pix_idx,
             l1c_L1qcOpts opts, l1c_LBResult *lb_res);

int l1c_breg_anistropic_TV(l1c_int n, l1c_int m, double *uk, double *f,
                       double mu, double tol, int max_iter, int max_jac_iter);


int l1c_dct2_setup(l1c_int Nx, l1c_int Mx, l1c_int Ny, l1c_int *pix_mask_idx,  l1c_AxFuns *ax_funs);

int l1c_dct1_setup(l1c_int Nx, l1c_int Ny, l1c_int *pix_mask_idx, l1c_AxFuns *ax_funs);

int l1c_setup_dct_transforms(l1c_int Nx, l1c_int Mx, l1c_int Ny,
                             l1c_int *pix_idx, l1c_AxFuns *ax_funs);

int l1c_setup_matrix_transforms(l1c_int n, l1c_int m, double *A, l1c_AxFuns *ax_funs);


#define L1C_INFEASIBLE_START  (1U << 1)
#define L1C_OUT_OF_MEMORY     (1U << 3)
#define L1C_DCT_INIT_FAILURE  (1U << 5)
#define L1C_FILE_READ_FAILURE (1U << 7)
#define L1C_CGSOLVE_FAILURE   (1U << 9)
#define L1C_INVALID_ARGUMENT  (1U << 11)


#endif
