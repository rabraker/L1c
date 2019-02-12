/* Dont compile if we dont have fftw */
#include "config.h"

#ifdef _USEFFTW3_

#include <stdlib.h>
#include <fftw3.h>
#include <omp.h>
#include <math.h>

#include "dct.h"
#include "l1c_common.h"


static fftw_plan dct_plan_MtEty;
static fftw_plan dct_plan_EMx;

static double * restrict dct_Ety_sparse; // Product of E^T*y, input to DCT.
static double * restrict dct_x;
static double * restrict dct_y;

static l1c_int *dct_pix_mask_idx;
static l1c_int dct_Ny;
static l1c_int dct_Nx;
static double dct_root_1_by_2N;

/**
   If return is not zero, do not call dct_destroy();
*/
int dct_setup(l1c_int Nx, l1c_int Ny, l1c_int *pix_mask_idx){
  int status = 0;
  fftw_init_threads();
  int n_proc = omp_get_num_procs();
  int n_thread = omp_get_max_threads();
  printf("num procs: %d, num thread: %d\n", n_proc, n_thread);
  fftw_plan_with_nthreads(n_thread);

  l1c_int i=0;
  dct_Ety_sparse = malloc_double(Nx);
  dct_x = malloc_double(Nx);
  dct_y = malloc_double(Nx);
  if (!dct_x || !dct_y || !dct_Ety_sparse){
    fprintf(stderr, "Error allocating memory in dct_setup");
    status = L1C_OUT_OF_MEMORY;
    goto fail;
  }

  for (i=0; i<Nx; i++){
    dct_Ety_sparse[i] = 0;
    dct_x[i] = 0;
    dct_y[i] = 0;
  }

  dct_Nx = Nx;
  dct_Ny = Ny;
  dct_pix_mask_idx = pix_mask_idx;

  dct_root_1_by_2N = sqrt(1.0 / ( (double) dct_Nx * 2)); // Normalization constant.

  // FFTW_PATIENT | FFTW_DESTROY_INPUT| FFTW_PRESERVE_INPUT; | FFTW_PRESERVE_INPUT;
  // FFTW_PATIENT seems to give us about 2.5 seconds for the 512x512 test image.
  unsigned flags = FFTW_PATIENT;

  fftw_r2r_kind dct_kind_MtEty = FFTW_REDFT10; // DCT-II, "the" DCT
  fftw_r2r_kind dct_kind_EMx   = FFTW_REDFT01; // DCT-III, “the” IDCT

  dct_plan_MtEty = fftw_plan_r2r_1d(Nx, dct_Ety_sparse, dct_x, dct_kind_MtEty, flags);
  dct_plan_EMx = fftw_plan_r2r_1d(Nx, dct_x, dct_y, dct_kind_EMx, flags);
  if ( !dct_plan_MtEty || !dct_plan_EMx){
    fprintf(stderr, "Failed to initialize FFTW3 dct plans.\n");
    status = L1C_DCT_INIT_FAILURE;
    goto fail;
  }

  return status;

 fail:
  dct_destroy();
  return status;
}

void dct_destroy(){
  fftw_free(dct_Ety_sparse);
  fftw_free(dct_x);
  fftw_free(dct_y);

  fftw_destroy_plan(dct_plan_EMx);
  fftw_destroy_plan(dct_plan_MtEty);
}


void dct_idct(double *x){
  // double *x_ = __builtin_assume_aligned(x, DALIGN);

  cblas_dcopy(dct_Nx, x, 1, dct_x, 1);
  dct_x[0] = dct_x[0] * sqrt(2.0); //proper scaling

  fftw_execute_r2r(dct_plan_EMx, dct_x, x);
  // normalize
  for (int i=0; i<dct_Nx; i++){
    x[i] =x[i] * dct_root_1_by_2N;
  }

}

void dct_EMx_new(double *x_fftw, double * restrict y){
  /* Compute y = E * M *dct_x, where M is the IDCT, and E is the subsampling matrix.

     --x_fftw should have size dct_Nx.
     -- x_fftw MUST have been allocated with fftw_alloc_real;
     -- y should have size at least dct_Ny.
     On exit, the first dct_Ny entries of y will contain the result of E * M *x.
  */
  l1c_int i=0;
  double *y_a = __builtin_assume_aligned(y, DALIGN);
  double *dcty_a = __builtin_assume_aligned(dct_y, DALIGN);

  // Leave the input vector unchanged.
  double x0_tmp = x_fftw[0];
  x_fftw[0] = x_fftw[0] * sqrt(2.0);

  fftw_execute_r2r(dct_plan_EMx, x_fftw, dcty_a);
  x_fftw[0] = x0_tmp;   /* Undo the scaling operation. */

  // Apply the subsampling operation. Assume that the indexes are ordered, and
  // do this in place
  for (i=0; i<dct_Ny; i++){
    y_a[i] = dcty_a[ dct_pix_mask_idx[i]] * dct_root_1_by_2N;
  }

}


void dct_MtEty( double * restrict y, double * restrict x){
  /*Apply x = M^T * E^T * y.
    The results are normalized to matlabs convention. That is, divided by 1/sqrt(2*N) and
    x[0] <- x[0] * /sqrt(2).

    -- y should have dimension at least dct_Ny
    -- x should have dimension at least dct_Nx
    -- neither x nor y are required to be allocated by fftw.
   */
  double *y_a = __builtin_assume_aligned(y, DALIGN);
  double *x_a = __builtin_assume_aligned(x, DALIGN);
  double *dctx_a = __builtin_assume_aligned(dct_x, DALIGN);
  double *Ety_a = __builtin_assume_aligned(dct_Ety_sparse, DALIGN);

  l1c_int i=0;
  // E^T * y --> y_sparse. y_sparse should be set to all zeros, and pix_mask does not change.
  // Returns a pointer to the array containing the result, which has length N.
  for (i=0; i<dct_Ny; i++){
    // dct_Ety_sparse[dct_pix_mask_idx[i]] = y[i] * dct_root_1_by_2N;
    Ety_a[dct_pix_mask_idx[i]] = y_a[i];
  }

  // Result contained in dct_x.
  fftw_execute(dct_plan_MtEty); //  M^T * dct_y_sparse

  dct_x[0] = dct_x[0]/sqrt(2.0);
  for(i=0; i<dct_Nx; i++){
    x_a[i] = dctx_a[i]  * dct_root_1_by_2N;
  }

}

void dct_MtEt_EMx_new(double * restrict x, double * restrict z){
  /* Performs the Multiplication z = (EM)^T * (EM) * x
     on the supplied vector x.

     -- x should have dimension at least dct_Nx.
     -- z should have dimension at least dct_Nx
   */

  // dct_EMx_new(x, dct_y); //writes output y to global, which belongs to fftw.
  // dct_MtEty(dct_y, z);

  /* It should be faster to implement E^T*Ex in shot
   */
  double *dcty_a = __builtin_assume_aligned(dct_y, DALIGN);
  double *dctEty_a = __builtin_assume_aligned(dct_Ety_sparse, DALIGN);
  double *dctx_a = __builtin_assume_aligned(dct_x, DALIGN);
  double *z_a = __builtin_assume_aligned(z, DALIGN);

  l1c_int i = 0;
  // double one_by_2N = 1.0 / ( (double) dct_Nx * 2);

  double x0_tmp = x[0];
  x[0] = x[0] * sqrt(2.0);
  fftw_execute_r2r(dct_plan_EMx, x, dct_y);
  x[0] = x0_tmp;

  /* Apply the subsampling operation. E^T*Ey. We could hold off on normalization:
  However, this seems to reduce accuracy, and gives us a larger residual (6097.209 vs 6096.907)
  yet takes 1000 more CG iterations in the end (for the 512 x 512 test image. Presumeably, the second DCT is less accurate and we we dont get as close to, e.g., dct(x) == M * x
  */

  for (i=0; i<dct_Ny; i++){
    dctEty_a[ dct_pix_mask_idx[i]] = dcty_a[ dct_pix_mask_idx[i]] * dct_root_1_by_2N;
  }

  // Result contained in dct_x.
  fftw_execute(dct_plan_MtEty); // This is the M^T * dct_y_sparse

  dct_x[0] = dct_x[0]/sqrt(2.0);
  for(i=0; i<dct_Nx; i++){
    z_a[i] = dctx_a[i]  * dct_root_1_by_2N;
  }

}

#else

typedef int make_iso_compilers_happy;

#endif
