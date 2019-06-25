/* Dont compile if we dont have fftw */
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <omp.h>
#include <math.h>

#include "cblas.h"
#include "l1c.h"

/*----- Forward declarations ----- */
static void dct_destroy();
static void dct_EMx(double *x_fftw, double *y);
static void dct_MtEty(double *y, double *x);
static void dct_MtEt_EMx(double *x_fftw, double *z);

static void dct_Mx(double *x, double *y);
static void dct_Mty(double *y, double *x);
static void dct_Ex(double *x, double *y);
static void dct_Ety(double *y, double *x);
static void Identity(double *x, double *z);

static fftw_plan dct_plan_MtEty;
static fftw_plan dct_plan_EMx;

static double * restrict dct_Ety_sparse; // Product of E^T*y, input to DCT.
static double * restrict dct_x;
static double * restrict dct_y;

static l1c_int *dct_pix_mask_idx;
static l1c_int dct_n;
static l1c_int dct_m;
static double dct_root_1_by_2N;


/**
 * @ingroup transforms Sub-sampled 1D DCT transform set.
 *
 * The function will populate an
 * l1c_AxFuns struct such that
 *
 * \f{align}{
 *    M &= \textrm{inverse discrete cosine transform}\\
 *    Ax &= EMx \\
 *    Aty &= M^TE^Ty \\
 *    AtAx &=  M^TE^TEMx)
 * \f}
 * where \f$M\f$ represents the inverse one dimensional
 * discrete cosine transform matrix and \f$E\f$ represents the subsampling.
 * The fields `Mt`, `E`, `Et` and `data` will be `NULL`.
 *
 * The dimension of the transform `A` is `n` by `m`.
 *
 *
 * The sub-sampling operator, \f$E\f$, is the identy matrix with rows removed.
 * In python notation:
 *
 * @code{python}
 *     E = I[pix_mask_idx, :]
 * @endcode
 *
 *
 * To de-allocate the memory reserved by dct1_setup(), call ax_funs.destroy().
 * Do not call ax_funs.destroy() if return value of dct1_setup() is non-zero.
 *
 * .. note:: Although this function treats the signal `x` as a 1D vector, it can also be
 *    used a 2D signal, e.g., image. In `c`, we hold a matrix in memory as a concatenated
 *    vector. L1C follows the convention of using row major order. Thus, an N by M
 *    grayscale image becomes an 1D N*M array, where the rows of the original matrix have
 *    been concatenated together. This is in contrast to the more standard (in linear
 *    algebra) notion of the \f$\textrm{vec}()\f$ operation, which concatenates a
 *    matrix in columm major order.
 *
 *
 * @param[in] n Number of elements in pix_mask_idx, and number of rows `A`.
 * @param[in] m Number of columns in the transform `A` and the number of
 *            of elements in the underlying signal.
 * @param[in]  pix_mask_idx indeces of locations of the subsampling.
 * @param[out] ax_funs A structure of function pointers which will be populated.
 *             On successfull exit, The fields Ax, Aty, AtAx, destroy, and M
 *             will be non-null. The fields MT, E, and ET will be null.
 *
 * @return  0 if succesfull. If unsuccesfull, returns L1C_DCT_INIT_FAILURE
 *          or L1C_OUT_OF_MEMORY.
 *
 *
 * @warning This function assumes that its inputs have already been sanitized. In
 *          particular, if `max(pix_mask_idx) > m`, then segfaults are likely to occur.
 */
int l1c_dct1_setup(l1c_int n, l1c_int m, l1c_int *pix_mask_idx, l1c_AxFuns *ax_funs){

  int status = 0;
#if defined(HAVE_FFTW3_THREADS)
  fftw_init_threads();
  // int n_proc = omp_get_num_procs();
  int n_thread = omp_get_max_threads();
  fftw_plan_with_nthreads(n_thread);
#endif

  l1c_int i=0;
  dct_Ety_sparse = l1c_malloc_double(m);
  dct_x = l1c_malloc_double(m);
  dct_y = l1c_malloc_double(m);
  if (!dct_x || !dct_y || !dct_Ety_sparse){
    fprintf(stderr, "Error allocating memory in dct_setup");
    status = L1C_OUT_OF_MEMORY;
    goto fail;
  }

  for (i=0; i<m; i++){
    dct_Ety_sparse[i] = 0;
    dct_x[i] = 0;
    dct_y[i] = 0;
  }

  dct_m = m;
  dct_n = n;
  dct_pix_mask_idx = pix_mask_idx;

  dct_root_1_by_2N = sqrt(1.0 / ( (double) dct_m * 2));

  /* Note: except for c2r and hc2r, the default is FFTW_PRESERVE_INPUT.
     It could be worth trying FFTW_DESTROY_INPUT, since this can potentially
     allow more efficient algorithms.

     FFTW_PATIENT seems to give us about 2.5 seconds for the 512x512 test image.
   */
  unsigned flags = FFTW_PATIENT;

  fftw_r2r_kind dct_kind_MtEty = FFTW_REDFT10; // DCT-II, "the" DCT
  fftw_r2r_kind dct_kind_EMx   = FFTW_REDFT01; // DCT-III, “the” IDCT

  dct_plan_MtEty = fftw_plan_r2r_1d(m, dct_Ety_sparse, dct_x, dct_kind_MtEty, flags);
  dct_plan_EMx = fftw_plan_r2r_1d(m, dct_x, dct_y, dct_kind_EMx, flags);
  if ( !dct_plan_MtEty || !dct_plan_EMx){
    fprintf(stderr, "Failed to initialize FFTW3 dct plans.\n");
    status = L1C_DCT_INIT_FAILURE;
    goto fail;
  }

  ax_funs->n = n;
  ax_funs->m = m;
  ax_funs->p = m;
  ax_funs->q = m;
  ax_funs->norm_W = 1.0;

  ax_funs->Ax = dct_EMx;
  ax_funs->Aty = dct_MtEty;
  ax_funs->AtAx = dct_MtEt_EMx;

  ax_funs->Mx = dct_Mx;
  ax_funs->Mty = dct_Mty;
  ax_funs->Rx = dct_Ex;
  ax_funs->Rty = dct_Ety;

  ax_funs->Wtx = Identity;
  ax_funs->Wz = Identity;

  ax_funs->destroy = dct_destroy;
  ax_funs->data = NULL;

  return status;

 fail:
  dct_destroy();
  return status;
}

static void dct_destroy(){
  l1c_free_double(dct_Ety_sparse);
  l1c_free_double(dct_x);
  l1c_free_double(dct_y);
  fftw_destroy_plan(dct_plan_EMx);
  fftw_destroy_plan(dct_plan_MtEty);
#if defined(HAVE_FFTW3_THREADS)
  fftw_cleanup_threads();
#endif
}

static void Identity(double *x, double *z){
  cblas_dcopy(dct_m, x, 1, z, 1);
}


static void dct_Mx(double *x, double *y){
  // double *x_ = __builtin_assume_aligned(x, DALIGN);

  x[0] = x[0] * sqrt(2.0); //proper scaling

  fftw_execute_r2r(dct_plan_EMx, x, y);
  // normalize
  for (int i=0; i<dct_m; i++){
    y[i] =y[i] * dct_root_1_by_2N;
  }

  x[0] = x[0] / sqrt(2.0); //undo scaling
}

static void dct_Mty( double *y, double *x){
  /*Apply x = M^T * y.
    The results are normalized to matlabs convention. That is, divided by 1/sqrt(2*N) and
    x[0] <- x[0] * /sqrt(2).

    -- y should have dimension at least dct_m.
       Should be aligned to 16 byte boundary for FFTW
    -- x should have dimension at least dct_m. Should be aligned
       to DALIGN boundary, for us.
   */
  double *x_a = __builtin_assume_aligned(x, DALIGN);

  l1c_int i=0;
  // Result contained in x_a.
  fftw_execute_r2r(dct_plan_MtEty, y, x_a); //  M^T * dct_y_sparse

  x_a[0] = x_a[0]/sqrt(2.0);
  for(i=0; i<dct_m; i++){
    x_a[i] = x_a[i]  * dct_root_1_by_2N;
  }

}


static void dct_EMx(double *x, double * restrict y){
  /* Compute y = E * M *dct_x, where M is the IDCT, and E is the subsampling matrix.
     --x should have size dct_m.
     -- x and y should have been allocated with l1c_malloc_double;
     -- y should have size at least dct_n.
     On exit, the first dct_n entries of y will contain the result of E * M *x.
  */
  l1c_int i=0;
  double *y_a = __builtin_assume_aligned(y, DALIGN);
  double *dcty_a = __builtin_assume_aligned(dct_y, DALIGN);

  // Leave the input vector unchanged.
  double x0_tmp = x[0];
  x[0] = x[0] * sqrt(2.0);

  fftw_execute_r2r(dct_plan_EMx, x, dcty_a);
  x[0] = x0_tmp;   /* Undo the scaling operation. */

  // Apply the subsampling operation. Assume that the indexes are ordered, and
  // do this in place
  for (i=0; i<dct_n; i++){
    y_a[i] = dcty_a[ dct_pix_mask_idx[i]] * dct_root_1_by_2N;
  }

}


static void dct_MtEty( double * restrict y, double * restrict x){
  /*Apply x = M^T * E^T * y.
    The results are normalized to matlabs convention. That is, divided by 1/sqrt(2*N) and
    x[0] <- x[0] * /sqrt(2).

    -- y should have dimension at least dct_n
    -- x should have dimension at least dct_m
    -- neither x nor y are required to be allocated by fftw.
   */
  double *y_a = __builtin_assume_aligned(y, DALIGN);
  double *x_a = __builtin_assume_aligned(x, DALIGN);
  double *dctx_a = __builtin_assume_aligned(dct_x, DALIGN);
  double *Ety_a = __builtin_assume_aligned(dct_Ety_sparse, DALIGN);

  l1c_int i=0;
  // E^T * y --> y_sparse. y_sparse should be set to all zeros, and pix_mask does not change.
  // Returns a pointer to the array containing the result, which has length N.
  for (i=0; i<dct_n; i++){
    // dct_Ety_sparse[dct_pix_mask_idx[i]] = y[i] * dct_root_1_by_2N;
    Ety_a[dct_pix_mask_idx[i]] = y_a[i];
  }

  // Result contained in dct_x.
  fftw_execute(dct_plan_MtEty); //  M^T * dct_y_sparse

  dct_x[0] = dct_x[0]/sqrt(2.0);
  for(i=0; i<dct_m; i++){
    x_a[i] = dctx_a[i]  * dct_root_1_by_2N;
  }

}

static void dct_MtEt_EMx(double * restrict x, double * restrict z){
  /* Performs the Multiplication z = (EM)^T * (EM) * x
     on the supplied vector x.

     -- x should have dimension at least dct_m.
     -- z should have dimension at least dct_m
   */

  // dct_EMx(x, dct_y); //writes output y to global, which belongs to fftw.
  // dct_MtEty(dct_y, z);

  /* It should be faster to implement E^T*Ex in shot
   */
  double *dcty_a = __builtin_assume_aligned(dct_y, DALIGN);
  double *dctEty_a = __builtin_assume_aligned(dct_Ety_sparse, DALIGN);
  double *dctx_a = __builtin_assume_aligned(dct_x, DALIGN);
  double *z_a = __builtin_assume_aligned(z, DALIGN);

  l1c_int i = 0;
  // double one_by_2N = 1.0 / ( (double) dct_m * 2);

  double x0_tmp = x[0];
  x[0] = x[0] * sqrt(2.0);
  fftw_execute_r2r(dct_plan_EMx, x, dct_y);
  x[0] = x0_tmp;

  /* Apply the subsampling operation. E^T*Ey. We could hold off on normalization:
  However, this seems to reduce accuracy, and gives us a larger residual (6097.209 vs 6096.907)
  yet takes 1000 more CG iterations in the end (for the 512 x 512 test image. Presumeably, the second DCT is less accurate and we we dont get as close to, e.g., dct(x) == M * x
  */

  for (i=0; i<dct_n; i++){
    dctEty_a[ dct_pix_mask_idx[i]] = dcty_a[ dct_pix_mask_idx[i]] * dct_root_1_by_2N;
  }

  // Result contained in dct_x.
  fftw_execute(dct_plan_MtEty); // This is the M^T * dct_y_sparse

  dct_x[0] = dct_x[0]/sqrt(2.0);
  for(i=0; i<dct_m; i++){
    z_a[i] = dctx_a[i]  * dct_root_1_by_2N;
  }

}


/* Apply the subsampling y = E * x,
   i.e., y = x[pix_idx]
*/
static void dct_Ex(double *x, double *y){
  for(int i=0; i<dct_n; i++){
    y[i] = x[dct_pix_mask_idx[i]];
  }
}


/*
  Apply the adjoint of the sub-sampling operation,
  x = E^T * y
  I.e., x = zeros(n), x[pix_idx] = y.

  x should already be initialized.
*/
static void dct_Ety(double *y, double *x){

  cblas_dscal(dct_m, 0, x, 1);

  for (int i=0; i<dct_n; i++){
    x[dct_pix_mask_idx[i]] = y[i];
  }
}
