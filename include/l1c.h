#ifndef _L1C_
#define _L1C_
/*
 * @file l1c.h
 *
 * The public API for l1c.
 */

/**
 * Integer type for libl1c. For now, this is int, but in the future
 * it may be convienent to change it.
 */
typedef int l1c_int;

/**
 * \defgroup memory Functions for dealing with memory allocation that are general to all
 * of l1c. Functions for managing memory, and in particular, aligned allocation needed
 * for vectorized math.
 * @{
 */

/** How many bytes l1c_malloc_double will use for alignment.
    The default is 64. This should be sufficient for up to AVX512. For SSE2
    and AVX2, you could redifine this 32.
 */
#ifndef DALIGN
#define DALIGN 64
#endif

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
void l1c_free_double(double* dptr);

void l1c_free_double_2D(int nrow, double** ddptr);
/**@}*/

/** @ingroup nesta
 * Instructs optimization to operate in synthesis mode, e.g.,
 * \f$\min_x ||x|| \text{ s.t. } ||RWx-b||<\sigma\f$

 * Instructs optimization to operate in analysis mode, e.g.,
 *  \f$\min_x ||W^Tx|| \text{ s.t. } ||Rx-b||<\sigma\f$
 */
typedef enum BpMode { analysis = (1U << 0), synthesis = (1U << 1) } BpMode;

/**
 * @defgroup transforms Functions for linear transformations.
 * @{
 */

typedef enum DctMode { dct1 = 1, dct2 = 2 } DctMode;

typedef struct _l1c_AxFuns l1c_AxFuns;

/** A struct of function pointers for linear transforms.
 *
 * Any new set of transforms must implement the base transforms, `.Ax`,
 * the adjoint, `.Aty`, and the normal transform `.AtAx`.
 *
 * These transforms can be set up for either basis pursuit synthesis or
 * analsysis problems. In either case, the problem can be described as
 *
 * \f$
 *    \min_{z} ||W^T x||_1 \text{s.t.} ~ ||Ax - b|| < \sigma
 * \f$
 *
 * --In synthesis mode, \f$ W = I\f$ and \f$ A = RM\f$.
 *
 * --In analysis mode, \f$ W\f$ is some dictionary (possibly overcomplete,
 *   or possibly the same as M), and \f$ A = R\f$.
 *
 *
 * For any specific set of transforms, if a one of these operators is not
 * needed, then the class should implement it as the identity operation.
 * The documentation for the class implementing l1cAxfuns should describe
 * which mode it is intended for.
 *
 * In general, we assume that
 * \f{align}{
 *   R & \in\mathbb{R}^{n\times m}\\
 *   M & \in\mathbb{R}^{m \times p}\\
 *   n \le & m \leq p
 * \f}
 *
 * @warning Upon entry, the output vector should be initialized.
 *
 * @see l1c_dct1_setup()
 *      l1c_dct2_setup()
 *      l1c_setup_dct_transforms()
 *      l1c_setup_matrix_transforms()
 *      l1c_setup_dctTV_transforms()
 *
 */
struct _l1c_AxFuns {
  /** Length of observation vector. n < m <=p*/
  l1c_int n;
  /** Length of true signal, number of columns in R. */
  l1c_int q;
  /** Number of columns in M and A. In analysis, q=m since M=I*/
  l1c_int m;
  /** Number of columns in W. p >=m. In synthesis, p=m, W=I*/
  l1c_int p;

  /** The spectral norm of the operator M (or at least an upper bound.
   Currently, only nesta requires this.*/
  double norm_W;

  /** Compute y=Ax */
  void (*Ax)(double* x, double* y);

  /** Compute x=Aty */
  void (*Aty)(double* y, double* x);

  /** Compute z=AtAx */
  void (*AtAx)(double* x, double* z);

  /** Compute x=Mx */
  void (*Mx)(double* x, double* y);

  /** Compute z=M^Tx */
  void (*Mty)(double* x, double* z);

  /** Compute x=Mx */
  void (*Wz)(double* z, double* x);

  /** Compute z=M^Tx */
  void (*Wtx)(double* x, double* z);

  /** Computes y = R * x */
  void (*Rx)(double* x, double* y);

  /** Computes x = R^T * x. x should already be
      initialized. */
  void (*Rty)(double* y, double* x);

  /**Release data allocated by the associated setup function.
     All implementations must define .`destroy`. */
  void (*destroy)(void);

  /** Reserved for future use. Hopefully, will be used
   * to make the optimizations re-entrant, rather than
   * file-global variables (as is currently done).*/
  void* data;
};

/** @}*/

/**
 * @defgroup lin_solve Routines for solving systems of linear equations.
 * @{*/
/** Struct containing artifacts of the cgsolve routine. */
typedef struct _l1c_CgResults {
  /** Residual */
  double cgres;
  /** Number of completed conjugate gradient iterations. */
  l1c_int cgiter;

} l1c_CgResults;

/**
 * Parameters for the conjugate gradient solver.
 */
typedef struct _l1c_CgParams {
  /** If 0, print nothing, if >0, print status every verbose-th iteration. */
  l1c_int verbose;
  /** Maximum number of solver iterations.*/
  l1c_int max_iter;
  /** Solver tolerance.*/
  double tol;
} l1c_CgParams;

/** @}*/

/**
 * @defgroup l1qc_lb l1-minimization with quadratic constraint.
 * @{*/

typedef struct _l1c_LBResult l1c_LBResult;
/**
 * Contains the results of an l1qc_newton() optimizations.
 */
struct _l1c_LBResult {
  /** Value of the objective function.*/
  double l1;
  /** Total number of newton interations, across all log-barrier iterations.*/
  int total_newton_iter;
  /** Total number of conjugate gradient interations, across all log-barrier
      and all Newton iterations.*/
  int total_cg_iter;
  /** Return status. 0 if the optimizations completed succesfully. Otherwise,
      see \ref l1c_errors */
  int status;
};

typedef struct _l1c_L1qcOpts l1c_L1qcOpts;
/**
 * Options which control the l1qc_dct() optimization.
 */
struct _l1c_L1qcOpts {
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
   * \f$ 0.5(\nabla f)^T \begin{bmatrix}d_x \\d_u\end{bmatrix} <
   * \texttt{newton}\_\texttt{tol} \f$ and the optimization continues into the next
   * log-barrier iteration.
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

  DctMode dct_mode;
};

/** @}*/

/**
 * @defgroup nesta Optimizations based on Nesterov's algorithm @cite
 * becker_nesta_2011
 * @{*/

/** Options for NESTA optimization.
 */
struct _l1c_NestaOpts {
  /** Smaller `mu` leads to better accuracy. Try 1e-5.*/
  double mu;
  /** Tolerance for termination criteria.*/
  double tol;
  /** Noise level of observations.*/
  double sigma;
  /**Number of continuation iterations. */
  int n_continue;
  /** Print nothing if 0.*/
  int verbose;
};

typedef struct _l1c_NestaOpts l1c_NestaOpts;

int l1c_nesta(l1c_int m,
              double* xk,
              l1c_int n,
              double* b,
              l1c_AxFuns ax_funs,
              l1c_NestaOpts opts);

/** @} */

int l1c_cgsolve(l1c_int N,
                double* x,
                double* b,
                double** Dwork,
                void (*AX_func)(l1c_int n, double* x, double* b, void* AX_data),
                void* AX_data,
                l1c_CgResults* cg_result,
                l1c_CgParams cg_params);

int l1c_cgsolve_diag_precond(
    l1c_int N,
    double* x,
    double* b,
    double* M_inv_diag,
    double** Dwork,
    void (*AX_func)(l1c_int n, double* x, double* b, void* AX_data),
    void* AX_data,
    l1c_CgResults* cg_result,
    l1c_CgParams cg_params);

l1c_LBResult l1c_l1qc_newton(l1c_int m,
                             double* x,
                             l1c_int n,
                             double* b,
                             l1c_L1qcOpts params,
                             l1c_AxFuns Ax_funs);

int l1qc_dct(int mrow,
             int mcol,
             double* x_out,
             int n,
             double* b,
             l1c_int* pix_idx,
             l1c_L1qcOpts opts,
             l1c_LBResult* lb_res);

int l1c_breg_anistropic_TV(l1c_int n,
                           l1c_int m,
                           double* uk,
                           double* f,
                           double mu,
                           double tol,
                           int max_iter,
                           int max_jac_iter);

int l1c_dct2_setup(
    l1c_int n, l1c_int mrow, l1c_int mcol, l1c_int* pix_mask_idx, l1c_AxFuns* ax_funs);

int l1c_dct1_setup(l1c_int n, l1c_int m, l1c_int* pix_mask_idx, l1c_AxFuns* ax_funs);

int l1c_setup_dct_transforms(l1c_int n,
                             l1c_int mrow,
                             l1c_int mcol,
                             DctMode dct_mode,
                             l1c_int* pix_idx,
                             l1c_AxFuns* ax_funs);

int l1c_setup_dctTV_transforms(l1c_int n,
                               l1c_int mrow,
                               l1c_int mcol,
                               double alp_v,
                               double alp_h,
                               DctMode dct_mode,
                               BpMode bp_mode,
                               l1c_int* pix_idx,
                               l1c_AxFuns* ax_funs);

int l1c_setup_matrix_transforms(l1c_int n, l1c_int m, double* A, l1c_AxFuns* ax_funs);

/**
 * @defgroup l1c_errors Error codes returned.
 * @{*/
/** Successfull return status*/
#define L1C_SUCCESS 0

/** Initial guess was infeasible. */
#define L1C_INFEASIBLE_START (1U << 1)

/** Call to malloc failed.*/
#define L1C_OUT_OF_MEMORY (1U << 3)

/** FFTW returned a null handle.*/
#define L1C_DCT_INIT_FAILURE (1U << 5)

/** Could not read file.*/
#define L1C_FILE_READ_FAILURE (1U << 7)

/** Conjugate gradient solver failed.*/
#define L1C_CGSOLVE_FAILURE (1U << 9)

/** E.g., N<0.*/
#define L1C_INVALID_ARGUMENT (1U << 11)

/** E.g., asking for analysis, but ax_funs.U==NULL*/
#define L1C_INCONSISTENT_ARGUMENTS (1U << 13)

/** @}*/

#endif
