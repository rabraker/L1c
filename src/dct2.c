#include "config.h"
#include <fftw3.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "cblas.h"
#include "l1c.h"

/* ----- Forward Declarations -----*/
static void dct2_destroy();
static void dct2_Mx(double* x, double* y);
static void dct2_Mty(double* y, double* x);
static void dct2_EMx(double* x_fftw, double* y);
static void dct2_MtEty(double* y, double* x);
static void dct2_MtEt_EMx(double* x_fftw, double* z);
static void dct2_Ex(double* x, double* y);
static void dct2_Ety(double* y, double* x);
static void Identity(double* y, double* x);

/*-------- Globals for dct2. ------------- */
static fftw_plan dct2_plan_MtEty;
static fftw_plan dct2_plan_EMx;

static double* restrict dct2_Ety_sparse; // Product of E^T*y, input to DCT.
static double* restrict dct2_x;
static double* restrict dct2_y;

static l1c_int* dct2_pix_mask_idx;
static l1c_int dct2_n;    // Number of measurements.
static l1c_int dct2_mrow; // Rows of image
static l1c_int dct2_mcol; // Cols of image

static double dct2_root_1_by_4NM;

/**
 * @ingroup transforms Sub-sampled 2D DCT transform set in synthesis mode.
 *
 * The function will populate an
 * l1c_AxFuns struct such that
 *
 * \f{align}{
 *    M &= \textrm{inverse discrete cosine transform}\\
 *    Mt &= \textrm{discrete cosine transform}\\
 *    Ax &= EM(x) \\
 *    Aty &= M^T(E^Ty) \\
 *    AtAx &=  M^T(E^TEM(x)) \\
 *    W=W^T&= I
 * \f}
 *
 * where \f$M(x)\f$ represents the inverse two dimensional
 * discrete cosine transform and \f$E\f$ represents the subsampling operation.
 * The fields `E`, `Et` and `data` will be `NULL`.
 *
 * The effective size of `A` is `n` by `mrow*mcol`.
 *
 * Recall that, for an `mrow` by `mcol` matrix \f$ X\f$, the 2D DCT is given by
 *
 * \f{equation}{
 *    \hat{X}=M_{mrow} X M_{mcol}
 * \f}
 *
 * On the other hand, the subsampling
 * operator operates on a vector. Specifically, \f$E\f$, is the identy matrix
 * with rows removed. In python notation:
 *
 * @code{python}
 *    E = I[pix_mask_idx,:]
 * @endcode
 *
 *
 * Conceptually, the operation
 *
 * \f$EM(X)\f$ can be thought of as
 *
 * \f{equation}{
 *    EM(x)= E\textrm{vec}\left(  (M_{mrow} X M_{mcol})^T  \right).
 * \f}
 *
 * The transpose is necessary because the `vec()` operator concatenates
 * a matrix in column major order, but, in L1C, we follow the standard
 * `C` convention of row major order. In reality, the code does not do any
 * of the re-shaping implied above, because the vector and matrix are represnted
 * the same way in memory. So an alternative description is that
 *
 * \f{equation}{n
 *    M(x)=M_{mcol}\otimes M^T_{mrow} \textrm{vec}(X^T)
 * \f}
 *
 * So that \f$M\f$ is an `(mrow * mcol)` by `(mrow * mcol)` matrix.
 *
 * To de-allocate the memory reserved by dct2_setup(), call ax_funs.destroy().
 * Do not call ax_funs.destroy() if return value of dct2_setup() is non-zero.
 *
 * @param[in] mrow Number of rows of the underlying signal.
 * @param[in] mcol Number of columns of the underlying signal.
 * @param[in]  n Number of elements in pix_mask_idx.
 * @param[in]  pix_mask_idx Indeces of locations of the subsampling. For both
 *              DCT1 and DCT2, this vector should be the same.
 * @param[out] ax_funs A structure of function pointers which will be populated.
 *                      On successfull exit, The fields Ax, Aty, AtAx, destroy,
 * and M will be non-null. The fields MT, E, and ET will be null.
 *
 *
 * @return   0 if succesfull. If unsuccesfull, returns L1C_DCT_INIT_FAILURE
 *           or L1C_OUT_OF_MEMORY.
 *
 * @warning This function assumes that its inputs have already been sanitized.
 * In particular, if `max(pix_mask_idx) > mrow*mcol`, then segfaults are likely
 * to occur.
 *
 * @note It should not be too hard to make this also work with analysis mode.
 *       However, that remains an outstanding goal.
 */
int l1c_dct2_setup(
    l1c_int n, l1c_int mrow, l1c_int mcol, l1c_int* pix_mask_idx, l1c_AxFuns* ax_funs) {
  int status = 0;
#if defined(HAVE_FFTW3_THREADS)
  fftw_init_threads();
  // int n_proc = omp_get_num_procs();
  int n_thread = omp_get_max_threads();
  fftw_plan_with_nthreads(n_thread);
#endif

  l1c_int i = 0;
  dct2_Ety_sparse = l1c_malloc_double(mrow * mcol);
  dct2_x = l1c_malloc_double(mrow * mcol);
  dct2_y = l1c_malloc_double(mrow * mcol);
  if (!dct2_x || !dct2_y || !dct2_Ety_sparse) {
    fprintf(stderr, "Error allocating memory in dct2_setup");
    status = L1C_OUT_OF_MEMORY;
    goto fail;
  }

  for (i = 0; i < mrow * mcol; i++) {
    dct2_Ety_sparse[i] = 0;
    dct2_x[i] = 0;
    dct2_y[i] = 0;
  }

  dct2_mrow = mrow;
  dct2_mcol = mcol;

  dct2_n = n;
  dct2_pix_mask_idx = pix_mask_idx;

  double den = ((double)dct2_mrow) * ((double)dct2_mcol) * 4;
  dct2_root_1_by_4NM = sqrt(1.0 / den); // Normalization constant.

  unsigned flags = FFTW_PATIENT;

  fftw_r2r_kind dct2_kind_MtEty = FFTW_REDFT10; // DCT-II, “the” DCT
  fftw_r2r_kind dct2_kind_EMx = FFTW_REDFT01;   // DCT-III, “the” IDCT

  dct2_plan_MtEty = fftw_plan_r2r_2d(
      mrow, mcol, dct2_Ety_sparse, dct2_x, dct2_kind_MtEty, dct2_kind_MtEty, flags);

  dct2_plan_EMx = fftw_plan_r2r_2d(mrow, mcol, dct2_x, dct2_y, dct2_kind_EMx, dct2_kind_EMx, flags);

  if (!dct2_plan_EMx || !dct2_plan_MtEty) {
    fprintf(stderr, "Failed to initialize FFTW3 dct plans.\n");
    status = L1C_DCT_INIT_FAILURE;
    goto fail;
  }

  ax_funs->n = n;
  ax_funs->m = mrow * mcol;
  ax_funs->p = mrow * mcol;
  ax_funs->q = mrow * mcol;
  ax_funs->norm_W = 1.0;
  ax_funs->Ax = dct2_EMx;
  ax_funs->Aty = dct2_MtEty;
  ax_funs->AtAx = dct2_MtEt_EMx;
  ax_funs->Mx = dct2_Mx;
  ax_funs->Mty = dct2_Mty;
  ax_funs->Rx = dct2_Ex;
  ax_funs->Rty = dct2_Ety;

  ax_funs->Wtx = Identity;
  ax_funs->Wz = Identity;

  ax_funs->destroy = dct2_destroy;
  ax_funs->data = NULL;

  return status;

fail:

  return status;
}

static void dct2_destroy() {
  l1c_free_double(dct2_Ety_sparse);
  l1c_free_double(dct2_x);
  l1c_free_double(dct2_y);

  fftw_destroy_plan(dct2_plan_EMx);
  fftw_destroy_plan(dct2_plan_MtEty);
#if defined(HAVE_FFTW3_THREADS)
  fftw_cleanup_threads();
#endif
}

static void Identity(double* y, double* x) { cblas_dcopy(dct2_mcol * dct2_mrow, y, 1, x, 1); }
/*
  Computes x = M^T * x
 */
static void dct2_Mty(double* y, double* x) {

  double one_by_root2 = 1.0 / sqrt(2);

  fftw_execute_r2r(dct2_plan_MtEty, y, x); //-->dct2_x

  // normalize by 1/sqrt(2*N) * 1/sqrt(2*M)
  cblas_dscal(dct2_mrow * dct2_mcol, dct2_root_1_by_4NM, x, 1);

  x[0] = x[0] * 0.5;
  // Across the columns of the first row.
  cblas_dscal(dct2_mcol - 1, one_by_root2, x + 1, 1);

  // Down the rows of the first column.
  cblas_dscal(dct2_mrow - 1, one_by_root2, x + dct2_mcol, dct2_mcol);
}

/*
  Computes y = M * x
 */
static void dct2_Mx(double* x, double* y) {

  double root2 = sqrt(2);
  cblas_dcopy(dct2_mrow * dct2_mcol, x, 1, dct2_x, 1);

  dct2_x[0] = dct2_x[0] * 2.0;

  // Across the columns of the first row.
  cblas_dscal(dct2_mcol - 1, root2, dct2_x + 1, 1);

  // Down the rows of the first column.
  cblas_dscal(dct2_mrow - 1, root2, dct2_x + dct2_mcol, dct2_mcol);

  fftw_execute_r2r(dct2_plan_EMx, dct2_x, y);

  // normalize by 1/sqrt(2*N) * 1/sqrt(2*M)
  cblas_dscal(dct2_mrow * dct2_mcol, dct2_root_1_by_4NM, y, 1);
}

static void dct2_EMx(double* x, double* y) {

  double root2 = sqrt(2);
  int i;
  cblas_dcopy(dct2_mrow * dct2_mcol, x, 1, dct2_x, 1);

  dct2_x[0] = dct2_x[0] * 2.0;

  // Across the columns of the first row.
  cblas_dscal(dct2_mcol - 1, root2, dct2_x + 1, 1);

  // Down the rows of the first column.
  cblas_dscal(dct2_mrow - 1, root2, dct2_x + dct2_mcol, dct2_mcol);

  fftw_execute_r2r(dct2_plan_EMx, dct2_x, dct2_y);

  // Apply the subsampling operation and
  // normalize by 1/sqrt(2*N) * 1/sqrt(2*M)
  for (i = 0; i < dct2_n; i++) {
    y[i] = dct2_y[dct2_pix_mask_idx[i]] * dct2_root_1_by_4NM;
  }
}

static void dct2_MtEty(double* y, double* x) {

  double one_by_root2 = 1.0 / sqrt(2);
  int i;
  for (i = 0; i < dct2_n; i++) {
    dct2_Ety_sparse[dct2_pix_mask_idx[i]] = y[i];
  }

  fftw_execute(dct2_plan_MtEty); //-->dct2_x

  // normalize by 1/sqrt(2*N) * 1/sqrt(2*M)
  for (i = 0; i < dct2_mrow * dct2_mcol; i++) {
    x[i] = dct2_x[i] * dct2_root_1_by_4NM;
  }

  x[0] = x[0] * 0.5;
  // Across the columns of the first row.
  cblas_dscal(dct2_mcol - 1, one_by_root2, x + 1, 1);

  // Down the rows of the first column.
  cblas_dscal(dct2_mrow - 1, one_by_root2, x + dct2_mcol, dct2_mcol);
}

static void dct2_MtEt_EMx(double* restrict x, double* restrict z) {
  (void)x;
  (void)z;

  double root2 = sqrt(2);
  double one_by_root2 = 1.0 / sqrt(2);
  int i;

  // /* ---------- EMx----------------------------*/
  cblas_dcopy(dct2_mrow * dct2_mcol, x, 1, dct2_x, 1);

  dct2_x[0] = dct2_x[0] * 2.0;

  // Across the columns of the first row.
  cblas_dscal(dct2_mcol - 1, root2, dct2_x + 1, 1);

  // Down the rows of the first column.
  cblas_dscal(dct2_mrow - 1, root2, dct2_x + dct2_mcol, dct2_mcol);

  fftw_execute_r2r(dct2_plan_EMx, dct2_x, dct2_y);

  // Apply the subsampling operation and
  // normalize by 1/sqrt(2*N) * 1/sqrt(2*M)
  for (i = 0; i < dct2_n; i++) {
    dct2_Ety_sparse[dct2_pix_mask_idx[i]] = dct2_y[dct2_pix_mask_idx[i]] * dct2_root_1_by_4NM;
    // dct2_Ety_sparse[dct2_pix_mask_idx[i]] = x[i];
  }

  /* ---------- MtEty----------------------------*/

  fftw_execute(dct2_plan_MtEty); //-->dct2_x

  // normalize by 1/sqrt(2*N) * 1/sqrt(2*M)
  for (i = 0; i < dct2_mrow * dct2_mcol; i++) {
    z[i] = dct2_x[i] * dct2_root_1_by_4NM;
  }

  z[0] = z[0] * 0.5;
  // Across the columns of the first row.
  cblas_dscal(dct2_mcol - 1, one_by_root2, z + 1, 1);

  // Down the rows of the first column.
  cblas_dscal(dct2_mrow - 1, one_by_root2, z + dct2_mcol, dct2_mcol);
}

/* Apply the subsampling y = E * x,
   i.e., y = x[pix_idx]
 */
static void dct2_Ex(double* x, double* y) {
  for (int i = 0; i < dct2_n; i++) {
    y[i] = x[dct2_pix_mask_idx[i]];
  }
}

/*
  Apply the adjoint of the sub-sampling operation,
  x = E^T * y
  I.e., x = zeros(n), x[pix_idx] = y.

  x should already be initialized.
 */
static void dct2_Ety(double* y, double* x) {

  cblas_dscal(dct2_mcol * dct2_mrow, 0, x, 1);

  for (int i = 0; i < dct2_n; i++) {
    x[dct2_pix_mask_idx[i]] = y[i];
  }
}
