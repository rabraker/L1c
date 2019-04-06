#ifndef __L1C__
#define __L1C__
/**
 * @file l1c.h
 *
 * The main API for l1c.
 */


/**
 * Integer type for libl1c. For now, this is int, but in the future
 * it may be convienent to change it.
 */
typedef int l1c_int;

/** How many bytes l1c_malloc_double will use for alignment.*/
#define DALIGN  64

/**
 *  Allocate an array of doubles aligned to a DALGIN byte boundary.
 *  @param[in] N Number of requested elements.
*/
double* l1c_malloc_double(int N);

/**
 * Free an array allocated by l1c_malloc_double
 */
void l1c_free_double(double *);

double** l1c_malloc_double_2D(l1c_int nrow, l1c_int ncol);

void l1c_free_double_2D(int nrow, double **dwork);


/**
 * A struct of function pointers for linear transforms.
 */
typedef struct l1c_AxFuns {
  /*Everybody must implement these.*/

  /** Compute y=Ax */
  void(*Ax)(double *x, double *y);
  /** Compute x=Aty */
  void(*Aty)(double *y, double *x);
  /** Compute z=AtAx */
  void(*AtAx)(double *x, double *z);

  /* Optional. Used by dct and dct2.*/
  /** Compute x=Mx*/
  void(*M)(double *x);
  /** Compute y=Mty*/
  void(*Mt)(double *y);

  /** Currently unused.*/
  void(*E)(double *x);
  /** Currently unused.*/
  void(*Et)(double *y);

  /**Release data allocated by setup.
     All implementations must define destroy.
  */
  void(*destroy)(void);

  /** Reserved for use.*/
  void *data;
}l1c_AxFuns;


/** Struct containing artifacts of the cgsolve routine. */
typedef struct l1c_CgResults{
  /** Residual */
  double cgres;
  /** Number of completed conjugate gradient iterations. */
  l1c_int cgiter;

} l1c_CgResults;


/**
 * Parameters for the conjugate gradient solver.
 */
typedef struct l1c_CgParams{
  /** If 0, print nothing, if >0, print status every verbose-th iteration. */
  l1c_int verbose;
  /** Maximum number of solver iterations.*/
  l1c_int max_iter;
  /** Solver tolerance.*/
  double tol;
} l1c_CgParams;


/**
 *  Contains the results of an l1qc_newton() optimizations.
 */
typedef struct l1c_LBResult{
  double l1;
  int    total_newton_iter;
  int    total_cg_iter;
  int    status;

}l1c_LBResult;


/**
 * Options which control the l1qc_dct() optimization.
 */
typedef struct l1c_L1qcOpts{
  /** The epsilon in \f$||Ax-b||_2 < \epsilon\f$  */
  double epsilon;
  /** log-barier parameter. */
  double mu;
  /** Log barrier tolerance. The number of log-barrier iterations
   * depends on this parameter.
   */
  double lbtol;
  /** This should be removed. */
  double tau;
  /** This should be removed. */
  int lbiter;
  /** The newton iterations stop if
   * \f$ 0.5\nabla f)^T \begin{bmatrix}d_x \\d_u\end{bmatrix} < newton_tol \f$
   * The optimization continues into the next log-barrier iteration.
  */
  double newton_tol;
  /** Maximum number of Newton iterations.*/
  int newton_max_iter;
  /** How much to print.*/
  int verbose;
  /** Stop Newton iterations when the relative difference in l1-norms
   * between iterations falls below l1_tol. The optimization will continue
   * into the next log-barrier iteration.
  */
  double l1_tol;
  /** See l1c_CgParams */
  double cg_tol;
  /** See l1c_CgParams */
  int cg_maxiter;
  /** See l1c_CgParams */
  int cg_verbose;
  /** Should the conjugate-gradient solver be warm-started?
   */
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


/** Initial guess was infeasible. */
#define L1C_INFEASIBLE_START  (1U << 1)

/** Call to malloc failed.*/
#define L1C_OUT_OF_MEMORY     (1U << 3)

/** FFTW returned a null handle.*/
#define L1C_DCT_INIT_FAILURE  (1U << 5)

/** Could not read file.*/
#define L1C_FILE_READ_FAILURE (1U << 7)

/** Conjugate gradient solver failed.*/
#define L1C_CGSOLVE_FAILURE   (1U << 9)

/** E.g., N<0.*/
#define L1C_INVALID_ARGUMENT  (1U << 11)


#endif
