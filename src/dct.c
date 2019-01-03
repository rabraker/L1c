#include <stdlib.h>
#include <fftw3.h>

#ifdef _USEMKL_
#include "fftw/fftw3_mkl.h"
#endif

#ifdef _USETHREADS_
// #include <omp.h>
// #include <fftw3_threads.h>
#include <pthread.h>
#endif

#include <math.h>
#include "dct.h"

#define PI  3.141592653589793
#include "l1qc_common.h"


static fftw_plan dct_plan_MtEty;
static fftw_plan dct_plan_EMx;

static double *dct_Ety_sparse; // Product of E^T*y, input to DCT.
static double *dct_x;
static double *dct_y;

static int *dct_pix_mask_idx;
static int dct_Ny;
static int dct_Nx;
static double dct_root_1_by_2N;


void dct_setup(int Nx, int Ny, int *pix_mask_idx){

  // local global.

  #ifdef _USETHREADS_
  fftw_init_threads();
  fftw_plan_with_nthreads(6);
  #endif

  int i=0;
  dct_Ety_sparse = fftw_alloc_real(Nx);
  dct_x = fftw_alloc_real(Nx);
  dct_y = fftw_alloc_real(Nx);
  for (i=0; i<Nx; i++){
    dct_Ety_sparse[i] = 0;
    dct_x[i] = 0;
    dct_y[i] = 0;
  }

  dct_Nx = Nx;
  dct_Ny = Ny;
  dct_pix_mask_idx = pix_mask_idx;

  dct_root_1_by_2N = sqrt(1.0 / ( (double) dct_Nx * 2)); // Normalization constant.

  // FFTW_PATIENT | FFTW_DESTROY_INPUT| FFTW_PRESERVE_INPUT;
  unsigned flags = FFTW_PATIENT;

  fftw_r2r_kind dct_kind_MtEty = FFTW_REDFT10; //computes an REDFT10 transform, a DCT-II
  fftw_r2r_kind dct_kind_EMx   = FFTW_REDFT01; //computes an REDFT01 transform, a DCT-III, “the” IDCT

  dct_plan_MtEty = fftw_plan_r2r_1d(Nx, dct_Ety_sparse, dct_x, dct_kind_MtEty, flags);
  dct_plan_EMx = fftw_plan_r2r_1d(Nx, dct_x, dct_y, dct_kind_EMx, flags);
}

void dct_destroy(){
  fftw_free(dct_Ety_sparse);
  fftw_free(dct_x);
  fftw_free(dct_y);

  fftw_destroy_plan(dct_plan_EMx);
  fftw_destroy_plan(dct_plan_MtEty);
}


double* idct_plain(double *x_fftw){
  /* Compute y = M *dct_x, where M is the IDCT, and E is the subsampling matrix.
     This performs the same function as dct_EMx, except that you can provide your own
     array.

     NOTE: x_fftw MUST have been allocated with fftw_alloc_real;

     On exit, the first N_pix_mask entries of y will contain the result of E * M *x.
  */
  /* Will fill y. Plan_AT has saved the pointer to x and y, so we
     dont have to supply, but should be updated by the caller.*/
  fftw_execute_r2r(dct_plan_EMx, x_fftw, dct_y);

  return dct_y;
}

void dct_EMx_new(double *x_fftw, double *y){
  /* Compute y = E * M *dct_x, where M is the IDCT, and E is the subsampling matrix.

     --x_fftw should have size dct_Nx.
     -- x_fftw MUST have been allocated with fftw_alloc_real;
     -- y should have size at least dct_Ny.
     On exit, the first dct_Ny entries of y will contain the result of E * M *x.
  */
  int i=0;

  // Leave the input vector unchanged.
  double x0_tmp = x_fftw[0];
  x_fftw[0] = x_fftw[0] * sqrt(2.0);

  fftw_execute_r2r(dct_plan_EMx, x_fftw, dct_y);

  // Apply the subsampling operation. Assume that the indexes are ordered, and
  // Do this in place
  for (i=0; i<dct_Ny; i++){
    y[i] = dct_y[ dct_pix_mask_idx[i]] * dct_root_1_by_2N;
  }
  /*Undo the scaling operation. */
  x_fftw[0] = x0_tmp;
}


void dct_MtEty( double *y, double *x){
  /*Apply x = M^T * E^T * y.
    The results are normalized to matlabs convention. That is, divided by 1/sqrt(2*N) and
    x[0] <- x[0] * /sqrt(2).

    -- y should have dimension at least dct_Ny
    -- x should have dimension at least dct_Nx
    -- neither x nor y are required to be allocated by fftw.
   */

  int i=0;
  // E^T * y --> y_sparse. y_sparse should be set to all zeros, and pix_mask does not change.
  // Returns a pointer to the array containing the result, which has length N.
  for (i=0; i<dct_Ny; i++){
    // dct_Ety_sparse[dct_pix_mask_idx[i]] = y[i] * dct_root_1_by_2N;
    dct_Ety_sparse[dct_pix_mask_idx[i]] = y[i];
  }
  fftw_execute(dct_plan_MtEty); // This is the M^T * dct_y_sparse

  // Result contained in dct_x. Return a pointer to that array. Normalize first coef to matlab convention
  dct_x[0] = dct_x[0]/sqrt(2.0);
  for(i=0; i<dct_Nx; i++){
    x[i] = dct_x[i]  * dct_root_1_by_2N;
  }

}

void dct_MtEt_EMx_new(double *x_fftw, double *z){
  /* Performs the Multiplication z = (EM)^T * (EM) * x
     on the supplied vector x.

     -- x should have dimension at least dct_Nx.
        NOTE: you MUST have allocated x with fftw_alloc_real()!
     -- z should have dimension at least dct_Nx
   */

  dct_EMx_new(x_fftw, dct_y); //writes output y to global, which belongs to fftw.
  dct_MtEty(dct_y, z);


}
