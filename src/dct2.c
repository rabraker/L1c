#include "config.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <omp.h>
#include <math.h>

#include "cblas.h"
#include "l1c.h"

/* ----- Forward Declarations -----*/
static void dct2_destroy();
static void dct2_idct(double *x_fftw);
static void dct2_dct(double *x);
static void dct2_EMx(double *x_fftw, double *y);
static void dct2_MtEty(double *y, double *x);
static void dct2_MtEt_EMx(double *x_fftw, double *z);

/*-------- Globals for dct2. ------------- */
static fftw_plan dct2_plan_MtEty;
static fftw_plan dct2_plan_EMx;

static double * restrict dct2_Ety_sparse; // Product of E^T*y, input to DCT.
static double * restrict dct2_x;
static double * restrict dct2_y;

static l1c_int *dct2_pix_mask_idx;
static l1c_int dct2_Ny;
static l1c_int dct2_N; // Rows
static l1c_int dct2_M; // Cols

static double dct2_root_1_by_4NM;


/**
 * The function will populate an
 * l1c_AxFuns struct such that
 *
 * \f{align}{
 *    M &= \textrm{inverse discrete cosine transform}\\
 *    Mt &= \textrm{discrete cosine transform}\\
 *    Ax &= EM(x) \\
 *    Aty &= M^T(E^Ty) \\
 *    AtAx &=  M^T(E^TEM(x))
 * \f}
 *
 * where \f$M(x)\f$ represents the inverse two dimensional
 * discrete cosine transform and \f$E\f$ represents the subsampling operation.
 * The fields `E`, `Et` and `data` will be `NULL`.
 *
 *
 * Recall that, for an n_x by m_x matrix X, the 2D DCT is given by
 *
 * \f{equation}{
 *    \hat{X}=M_{n_x} X M_{m_x}
 * \f}
 *
 * where \f$ X\in\mathbb{R}^{n_x \times m_x}\f$. On the other hand, the subsampling
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
 *    EM(x)= E\textrm{vec}\left(  (M_{n_x} X M_{m_x})^T  \right).
 * \f}
 *
 * The transpose is necessary because the `vec()` operator concatenates
 * a matrix in column major order, but, in L1C, we follow the standard
 * `C` convention of row major order. In reality, the code does not do any
 * of the re-shaping implied above, because the vector and matrix are represnted
 * the same way in memory. So an alternative description is that
 *
 * \f{equation}{
 *    M(x)=M_{m_x}\otimes M^T_{n_x} \textrm{vec}(X^T)
 * \f}
 *
 *
 * To de-allocate the memory reserved by dct2_setup(), call ax_funs.destroy().
 * Do not call ax_funs.destroy() if return value of dct2_setup() is non-zero.
 *
 * @param[in] Nx Number of rows of the underlying signal.
 * @param[in] Mx Number of columns of the underlying signal.
 *               To treat, e.g., an image as a 1D vectorized signal
 *               set Mx=1 and Nx = number_of_rows * number_of_columns.
 * @param[in]  Ny Number of elements in pix_mask_idx.
 * @param[in]  pix_mask_idx Indeces of locations of the subsampling. For both
 *              DCT1 and DCT2, this vector should be the same.
 * @param[out] ax_funs A structure of function pointers which will be populated.
 *                      On successfull exit, The fields Ax, Aty, AtAx, destroy, and M
 *                      will be non-null. The fields MT, E, and ET will be null.
 *
 *
 * @return   0 if succesfull. If unsuccesfull, returns L1C_DCT_INIT_FAILURE
 *           or L1C_OUT_OF_MEMORY.
 *
 * @warning This function assumes that its inputs have already been sanitized. In
 *          particular, if `max(pix_mask_idx) > Nx*Mx`, then segfaults are likely to occur.
 */
int l1c_dct2_setup(l1c_int Nx, l1c_int Mx, l1c_int Ny, l1c_int *pix_mask_idx,  l1c_AxFuns *ax_funs){
  int status = 0;
#if defined(HAVE_FFTW3_THREADS)
  fftw_init_threads();
  // int n_proc = omp_get_num_procs();
  int n_thread = omp_get_max_threads();
  fftw_plan_with_nthreads(n_thread);
#endif

  l1c_int i=0;
  dct2_Ety_sparse = l1c_malloc_double(Nx*Mx);
  dct2_x = l1c_malloc_double(Nx*Mx);
  dct2_y = l1c_malloc_double(Nx*Mx);
  if (!dct2_x || !dct2_y || !dct2_Ety_sparse){
    fprintf(stderr, "Error allocating memory in dct2_setup");
    status = L1C_OUT_OF_MEMORY;
    goto fail;
  }

  for (i=0; i<Nx*Mx; i++){
    dct2_Ety_sparse[i] = 0;
    dct2_x[i] = 0;
    dct2_y[i] = 0;
  }

  dct2_N = Nx;
  dct2_M = Mx;

  dct2_Ny = Ny;
  dct2_pix_mask_idx = pix_mask_idx;

  double den = ((double) dct2_N) * ((double)dct2_M) * 4 ;
  dct2_root_1_by_4NM = sqrt(1.0 / den); // Normalization constant.


  unsigned flags = FFTW_PATIENT;

  fftw_r2r_kind dct2_kind_MtEty = FFTW_REDFT10; // DCT-II, “the” DCT
  fftw_r2r_kind dct2_kind_EMx = FFTW_REDFT01; // DCT-III, “the” IDCT

  dct2_plan_MtEty = fftw_plan_r2r_2d(Nx, Mx, dct2_Ety_sparse, dct2_x,
                             dct2_kind_MtEty, dct2_kind_MtEty, flags);

  dct2_plan_EMx = fftw_plan_r2r_2d(Nx, Mx, dct2_x, dct2_y,
                                     dct2_kind_EMx, dct2_kind_EMx, flags);

  if ( !dct2_plan_EMx || !dct2_plan_MtEty){
    fprintf(stderr, "Failed to initialize FFTW3 dct plans.\n");
    status = L1C_DCT_INIT_FAILURE;
    goto fail;
  }

  ax_funs->Ax = dct2_EMx;
  ax_funs->Aty = dct2_MtEty;
  ax_funs->AtAx = dct2_MtEt_EMx;
  ax_funs->M = dct2_idct;
  ax_funs->Mt = dct2_dct;
  ax_funs->E = NULL;
  ax_funs->Et = NULL;

  ax_funs->destroy = dct2_destroy;
  ax_funs->data = NULL;

  return status;

 fail:

  return status;
}


static void dct2_destroy(){
  l1c_free_double(dct2_Ety_sparse);
  l1c_free_double(dct2_x);
  l1c_free_double(dct2_y);

  fftw_destroy_plan(dct2_plan_EMx);
  fftw_destroy_plan(dct2_plan_MtEty);
#if defined(HAVE_FFTW3_THREADS)
  fftw_cleanup_threads();
#endif
}


static void dct2_dct(double *x){

  double one_by_root2 = 1.0/sqrt(2);
  int i;
  cblas_dcopy(dct2_N * dct2_M, x, 1, dct2_Ety_sparse, 1);


  fftw_execute(dct2_plan_MtEty); //-->dct2_x

  // normalize by 1/sqrt(2*N) * 1/sqrt(2*M)
  for (i=0; i<dct2_N * dct2_M; i++){
    x[i] = dct2_x[i] * dct2_root_1_by_4NM;
  }

  x[0] = x[0]*0.5;
  // Across the columns of the first row.
  cblas_dscal(dct2_M-1, one_by_root2, x+1, 1);

  // Down the rows of the first column.
  cblas_dscal(dct2_N-1, one_by_root2, x+dct2_M, dct2_M);

}


static void dct2_idct(double *x){

  double root2 = sqrt(2);
  cblas_dcopy(dct2_N * dct2_M, x, 1, dct2_x, 1);


  dct2_x[0] = dct2_x[0]*2.0;

  // Across the columns of the first row.
  cblas_dscal(dct2_N-1, root2, dct2_x+1, 1);

  // Down the rows of the first column.
  cblas_dscal(dct2_N-1, root2, dct2_x+dct2_M, dct2_M);

  fftw_execute_r2r(dct2_plan_EMx, dct2_x, x);

  // normalize by 1/sqrt(2*N) * 1/sqrt(2*M)
  cblas_dscal(dct2_N*dct2_M, dct2_root_1_by_4NM, x, 1);

}


static void dct2_EMx(double *x, double *y){

  double root2 = sqrt(2);
  int i;
  cblas_dcopy(dct2_N * dct2_M, x, 1, dct2_x, 1);


  dct2_x[0] = dct2_x[0]*2.0;

  // Across the columns of the first row.
  cblas_dscal(dct2_M-1, root2, dct2_x+1, 1);

  // Down the rows of the first column.
  cblas_dscal(dct2_N-1, root2, dct2_x+dct2_M, dct2_M);

  fftw_execute_r2r(dct2_plan_EMx, dct2_x, dct2_y);

  // Apply the subsampling operation and
  // normalize by 1/sqrt(2*N) * 1/sqrt(2*M)
  for(i=0; i<dct2_Ny; i++){
    y[i] = dct2_y[dct2_pix_mask_idx[i]] * dct2_root_1_by_4NM;
  }

}


static void dct2_MtEty(double *y, double *x){

  double one_by_root2 = 1.0/sqrt(2);
  int i;
  for (i=0; i<dct2_Ny; i++){
    dct2_Ety_sparse[dct2_pix_mask_idx[i]] = y[i];
  }


  fftw_execute(dct2_plan_MtEty); //-->dct2_x

  // normalize by 1/sqrt(2*N) * 1/sqrt(2*M)
  for (i=0; i<dct2_N * dct2_M; i++){
    x[i] = dct2_x[i] * dct2_root_1_by_4NM;
  }

  x[0] = x[0]*0.5;
  // Across the columns of the first row.
  cblas_dscal(dct2_M-1, one_by_root2, x+1, 1);

  // Down the rows of the first column.
  cblas_dscal(dct2_N-1, one_by_root2, x+dct2_M, dct2_M);

}

static void dct2_MtEt_EMx(double * restrict x, double * restrict z){
  (void) x;
  (void) z;

  double root2 = sqrt(2);
  double one_by_root2 = 1.0/sqrt(2);
  int i;

  // /* ---------- EMx----------------------------*/
  cblas_dcopy(dct2_N * dct2_M, x, 1, dct2_x, 1);


  dct2_x[0] = dct2_x[0]*2.0;

  // Across the columns of the first row.
  cblas_dscal(dct2_M-1, root2, dct2_x+1, 1);

  // Down the rows of the first column.
  cblas_dscal(dct2_N-1, root2, dct2_x+dct2_M, dct2_M);

  fftw_execute_r2r(dct2_plan_EMx, dct2_x, dct2_y);

  // Apply the subsampling operation and
  // normalize by 1/sqrt(2*N) * 1/sqrt(2*M)
  for(i=0; i<dct2_Ny; i++){
    dct2_Ety_sparse[dct2_pix_mask_idx[i]] = dct2_y[dct2_pix_mask_idx[i]] * dct2_root_1_by_4NM;
    // dct2_Ety_sparse[dct2_pix_mask_idx[i]] = x[i];
  }

  /* ---------- MtEty----------------------------*/

  fftw_execute(dct2_plan_MtEty); //-->dct2_x

  // normalize by 1/sqrt(2*N) * 1/sqrt(2*M)
  for (i=0; i<dct2_N * dct2_M; i++){
    z[i] = dct2_x[i] * dct2_root_1_by_4NM;
  }

  z[0] = z[0]*0.5;
  // Across the columns of the first row.
  cblas_dscal(dct2_M-1, one_by_root2, z+1, 1);

  // Down the rows of the first column.
  cblas_dscal(dct2_N-1, one_by_root2, z+dct2_M, dct2_M);
}
