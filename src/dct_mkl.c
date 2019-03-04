#include "config.h"
#ifdef _USEMKL_ //disable the entire file if no mkl.

#include <stdlib.h>
#include <math.h>

#include "mkl_service.h"
#include "mkl_trig_transforms.h"
#include "l1c_common.h"

DFTI_DESCRIPTOR_HANDLE mkl_dct_handle=0;


static double *dct_tmp_x=NULL;
static double *dct_tmp_y=NULL;

static MKL_INT dct_Ny;
static MKL_INT dct_Nx;

static MKL_INT dct_ipar[128];
static double *dct_dpar=NULL;
static MKL_INT dct_stat;

static l1c_int *dct_pix_mask_idx;
static double dct_root_1_by_2N;
/*
  Standard BLAS interfaces (e.g., NetLib, ATLAS) use a standard int (32).

  OpenBlas ~can~ use long int (64) if OPENBLAS_USE64BITINT is defined.

  But for compatibility with as many as possible, we for now use a standard int.

  The int size in MKL depends on if one chooses ILP64 or LP64.
 */


/**
   If return is not zero, do not call dct_destroy();
 */
int dct_setup(l1c_int Nx, l1c_int Ny, l1c_int *pix_mask_idx){

  MKL_INT tt_type = MKL_STAGGERED_COSINE_TRANSFORM;
  l1c_int i=0;
  int status = 0;
  dct_Nx = (MKL_INT)Nx;
  dct_Ny = (MKL_INT)Ny;
  dct_root_1_by_2N = sqrt(1.0 / ( (double) dct_Nx * 2)); // Normalization constant.
  dct_pix_mask_idx = pix_mask_idx;

  dct_tmp_x = malloc_double(Nx+1);
  dct_tmp_y = malloc_double(Nx+1);
  dct_dpar = malloc_double(5*Nx/2 + 2);

  if (!dct_tmp_x || !dct_tmp_y || !dct_dpar){
    fprintf(stderr, "Error allocating memory in dct_setup\n");
    status = L1C_OUT_OF_MEMORY;
    goto fail;
  }

  for (i=0; i<Nx+1; i++){
    dct_tmp_x[i] = 0;
    dct_tmp_y[i] = 0;
  }


  d_init_trig_transform(&dct_Nx, &tt_type, dct_ipar, dct_dpar, &dct_stat);
  if (dct_stat){
    status = L1C_DCT_INIT_FAILURE;
    goto fail;
  }
  d_commit_trig_transform(dct_tmp_x, &mkl_dct_handle, dct_ipar, dct_dpar, &dct_stat);
  if (dct_stat){
    status = L1C_DCT_INIT_FAILURE;
    goto fail;
  }

  return status;

 fail:
  free_double(dct_tmp_x);
  free_double(dct_tmp_y);
  free_double(dct_dpar);
  return status;
}

void dct_destroy(){
  free_trig_transform(&mkl_dct_handle, dct_ipar, &dct_stat);

  free_double(dct_tmp_x);
  free_double(dct_tmp_y);
  free_double(dct_dpar);
}

void dct_idct(double * restrict x){
  double *x_ = __builtin_assume_aligned(x, DALIGN);

  x_[0] = x_[0] * sqrt(2.0); //proper scaling
  //dct_ipar[7] = 0;
  dct_ipar[10] = 1; // Don't normalize

  d_forward_trig_transform(x_, &mkl_dct_handle,
                            dct_ipar, dct_dpar, &dct_stat);
  // normalize
  for (int i=0; i<dct_Nx; i++){
    x_[i] =x_[i] * dct_root_1_by_2N;
  }

}

void dct_EMx(double * restrict x, double * restrict y){
  /* Compute y = E * M *dct_x, where M is the IDCT, and E is the subsampling matrix.

     -- x should have size dct_Nx.
     -- y should have size at least dct_Ny.
     On exit, the first dct_Ny entries of y will contain the result of E * M *x.
  */

  double *x_a = __builtin_assume_aligned(x, DALIGN);
  double *y_a = __builtin_assume_aligned(y, DALIGN);
  double *tmp_x = __builtin_assume_aligned(dct_tmp_x, DALIGN);

  l1c_int i=0;


  cblas_dcopy(dct_Nx, x_a, 1, tmp_x, 1);
  tmp_x[0] = tmp_x[0] * sqrt(2.0);  //proper scaling, so transforms are transposes.
  tmp_x[dct_Nx] = 0.0;           //Strange mkl requirement. dct_tmp_x should have size N+1


  dct_ipar[10] = 1; //Dont normalize the result, because the do it wrong.
  d_forward_trig_transform(tmp_x, &mkl_dct_handle, dct_ipar, dct_dpar, &dct_stat);

  //Apply the subsampling operation.
  for (i=0; i<dct_Ny; i++){
    y_a[i] = tmp_x[ dct_pix_mask_idx[i]] * dct_root_1_by_2N;
  }

}


void dct_MtEty( double * restrict y, double * restrict x){
  /*Apply x = M^T * E^T * y.
    -- y should have dimension at least dct_Ny
    -- x should have dimension at least dct_Nx
  */

  l1c_int i=0;
  double *x_a = __builtin_assume_aligned(x, DALIGN);
  double *y_a = __builtin_assume_aligned(y, DALIGN);
  double *tmp_x_a = __builtin_assume_aligned(dct_tmp_x, DALIGN);

  // E^T * y --> y_sparse. y_sparse should be set to all zeros, and pix_mask does not change.
  // Returns a pointer to the array containing the result, which has length N.

  for(i=0; i<dct_Nx+1; i++){
    tmp_x_a[i] = 0.0;
  }
  for(i=0; i<dct_Ny; i++){
    tmp_x_a[dct_pix_mask_idx[i]] = y_a[i];
  }


  dct_ipar[10] = 1; //Dont normalize
  d_backward_trig_transform(tmp_x_a, &mkl_dct_handle, dct_ipar, dct_dpar, &dct_stat);

  tmp_x_a[0] = tmp_x_a[0]/sqrt(2.0);

  cblas_daxpby(dct_Nx, dct_root_1_by_2N, tmp_x_a, 1, 0.0, x_a, 1);
}


void dct_MtEt_EMx(double *x_fftw, double *z){
  /* Performs the Multiplication z = (EM)^T * (EM) * x
     on the supplied vector x.

     -- x should have dimension at least dct_Nx.
     -- z should have dimension at least dct_Nx
  */

  dct_EMx(x_fftw, dct_tmp_y); //writes output y to global, which belongs to fftw.
  dct_MtEty(dct_tmp_y, z);


}

#else

typedef int make_iso_compilers_happy;

#endif
