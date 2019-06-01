#ifndef __L1C__
#define __L1C__
/*
 * @file l1c.h
 *
 * The main API for l1c.
 */


/**
 * Integer type for libl1c. For now, this is int, but in the future
 * it may be convienent to change it.
 */
typedef int l1c_int;


/**
 * \defgroup memory Functions for dealing with memory allocation that are general to all of l1c.
 * Functions for managing memory, and in particular, aligned allocation
 * needed for vectorized math.
 * @{
 */

/** How many bytes l1c_malloc_double will use for alignment.
    The default is 64. This should be sufficient for up to AVX512. For SSE2
    and AVX2, you could redifine this 32.
 */
#define DALIGN  64

/**
 * Allocate an array of doubles aligned to a DALGIN byte boundary.
 *
 * @param[in] N Number of requested elements. The returned pointer should be
 * freed with l1c_free_double().
 *
 * @return pointer to an array of doubles of length N.
 *
 * @see{l1c_calloc_double()}
*/
double* l1c_malloc_double(int N);

double* l1c_calloc_double(int N);

double** l1c_calloc_double_2D(l1c_int nrow, l1c_int ncol);

double** l1c_malloc_double_2D(l1c_int nrow, l1c_int ncol);


/**
 * Free an array allocated by l1c_malloc_double.
 *
 * @param[in] dptr The pointer to be freed.
 */
void l1c_free_double(double *dptr);

void l1c_free_double_2D(int nrow, double **ddptr);
/**@}*/



/**
 * @defgroup transforms Functions for linear transformations.
 * @{
 */

/** A struct of function pointers for linear transforms.
 *
 * Any new set of transforms must implement the base transforms, `.Ax`,
 * the adjoint, `.Aty`, and the normal transform `.AtAx`.
 *
 * Upon entry, the output vector should be initialized to any finite value.
 * This allows more efficient code in the matrix_transforms class, which is
 * based on cblas_dgemv, which does overwrite the output but rather updates
 * it in place. For example,
 *
 * @code
 * y = Ax + beta*y
 * @endcode
 *
 * The remaining function pointers are, at this point, optional, and not required
 * by the low level optimizations. However, often in CS, the transform `A` can be
 * decomposed as `A=EM`, where `M` is an othogonal, `m` by `m` matrix (say, a
 * DCT transform) and `E` is an `n` by `m` matrix, e.g., an identity with rows removed.
 * The remaining fields are there to provide access to those individual components.
 *
 * Currently, the `.data` field is unused, but at some point, it may be good to make the
 * optimizations re-entrant. This field is there as a placeholder so a transform class
 * can store its state in `.data`, rather than file-global variables (as is currently done).
 * @see l1c_dct1_setup()
 *      l1c_dct2_setup()
 *      l1c_setup_dct_transforms()
 *      l1c_setup_matrix_transforms()
 *
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
  void(*Mx)(double *x, double *y);
  /** Compute y=Mty*/
  void(*Mty)(double *y, double *x);

  /** Computes y = E * x
   */
  void(*Ex)(double *x, double *y);

  /** Computes x = E^T * x. x should already be
      initialized.
   */
  void(*Ety)(double *y, double *x);

  /**Release data allocated by the associated setup function.
     All implementations must define .`destroy`.
  */
  void(*destroy)(void);

  /** Reserved for future use.*/
  void *data;
}l1c_AxFuns;




/** @}*/

/**
 * @defgroup lin_solve Routines for solving systems of linear equations.
 * @{*/
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

/** @}*/



/**
 *  Contains the results of an l1qc_newton() optimizations.
 */
typedef struct l1c_LBResult{
  /** Value of the objective function.*/
  double l1;
  /** Total number of newton interations, across all log-barrier iterations.*/
  int    total_newton_iter;
  /** Total number of conjugate gradient interations, across all log-barrier
      and all Newton iterations.*/
  int    total_cg_iter;
  /** Return status. 0 if the optimizations completed succesfully. Otherwise,
      see \ref l1c_errors */
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

l1c_LBResult l1c_l1qc_newton(l1c_int m, double *x, l1c_int n, double *b,
                             l1c_L1qcOpts params, l1c_AxFuns Ax_funs);

int l1qc_dct(int mrow, int mcol, double *x_out, int n, double *b, l1c_int *pix_idx,
             l1c_L1qcOpts opts, l1c_LBResult *lb_res);

int l1c_breg_anistropic_TV(l1c_int n, l1c_int m, double *uk, double *f,
                       double mu, double tol, int max_iter, int max_jac_iter);


int l1c_dct2_setup(l1c_int mrow, l1c_int mcol, l1c_int n, l1c_int *pix_mask_idx,  l1c_AxFuns *ax_funs);

int l1c_dct1_setup(l1c_int m, l1c_int n, l1c_int *pix_mask_idx, l1c_AxFuns *ax_funs);


int l1c_setup_dct_transforms(l1c_int mrow, l1c_int mcol, l1c_int n,
                             l1c_int *pix_idx, l1c_AxFuns *ax_funs);

int l1c_setup_matrix_transforms(l1c_int n, l1c_int m, double *A, l1c_AxFuns *ax_funs);


/**
 * @defgroup l1c_errors Error codes returned.
 * @{*/
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

/** @}*/

#endif
