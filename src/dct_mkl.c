#ifdef _USEMKL_ //disable the entire file if no mkl.

#include <stdlib.h>
#include <math.h>

#include "mkl_service.h"
#include "mkl_trig_transforms.h"
#include "l1c_common.h"

DFTI_DESCRIPTOR_HANDLE mkl_dct_handle=0;


static double *dctmkl_tmp_x=NULL;
static double *dctmkl_tmp_y=NULL;

static l1c_int  dctmkl_Ny;
static l1c_int dctmkl_Nx;

static l1c_int dctmkl_ipar[128];
static double *dctmkl_dpar=NULL;
static l1c_int dctmkl_stat;

static l1c_int *dctmkl_pix_mask_idx;
static double dctmkl_root_1_by_2N;

int dctmkl_setup(l1c_int Nx, l1c_int Ny, l1c_int *pix_mask_idx){

  l1c_int tt_type = MKL_STAGGERED_COSINE_TRANSFORM;
  l1c_int i=0;

  dctmkl_Nx = Nx;
  dctmkl_Ny = Ny;
  dctmkl_root_1_by_2N = sqrt(1.0 / ( (double) dctmkl_Nx * 2)); // Normalization constant.
  dctmkl_pix_mask_idx = pix_mask_idx;

  dctmkl_tmp_x = malloc_double(Nx+1);
  dctmkl_tmp_y = malloc_double(Nx+1);
  dctmkl_dpar = malloc_double(5*Nx/2 + 2);

  if (!dctmkl_tmp_x | !dctmkl_tmp_y | !dctmkl_dpar){
    perror("Error allocating memory in dctmkl_setup");
    free_double(dctmkl_tmp_x);
    free_double(dctmkl_tmp_y);
    free_double(dctmkl_dpar);
    return 1;
  }

  for (i=0; i<Nx+1; i++){
    dctmkl_tmp_x[i] = 0;
    dctmkl_tmp_y[i] = 0;
  }


  d_init_trig_transform(&dctmkl_Nx, &tt_type, dctmkl_ipar, dctmkl_dpar, &dctmkl_stat);
  d_commit_trig_transform(dctmkl_tmp_x, &mkl_dct_handle, dctmkl_ipar, dctmkl_dpar, &dctmkl_stat);

  return 0;
}

void dctmkl_destroy(){
  free_trig_transform(&mkl_dct_handle, dctmkl_ipar, &dctmkl_stat);

  free_double(dctmkl_tmp_x);
  free_double(dctmkl_tmp_y);
  free_double(dctmkl_dpar);
}

void dctmkl_idct(double * restrict x, double * restrict y){
  double *x_ = __builtin_assume_aligned(x, 64);
  double *y_ = __builtin_assume_aligned(y, 64);

  for (l1c_int i=0; i<dctmkl_Nx; i++){
      y_[i] = x_[i];
    }
  // cblas_dcopy(dctmkl_Nx, x, 1, y, 1);
  dctmkl_ipar[7] = 0;
  dctmkl_ipar[10] = 1;

  d_forward_trig_transform(y_, &mkl_dct_handle,
                            dctmkl_ipar, dctmkl_dpar, &dctmkl_stat);
}

void dctmkl_EMx_new(double * restrict x, double * restrict y){
  /* Compute y = E * M *dctmkl_x, where M is the IDCT, and E is the subsampling matrix.

     --x_fftw should have size dctmkl_Nx.
     -- y should have size at least dctmkl_Ny.
     On exit, the first dctmkl_Ny entries of y will contain the result of E * M *x.
  */

  double *x_a = __builtin_assume_aligned(x, 64);
  double *y_a = __builtin_assume_aligned(y, 64);
  double *tmp_x = __builtin_assume_aligned(dctmkl_tmp_x, 64);
  // l1c_int idx_ = __builtin_assume_aligned(dctmkl_pix_mask_idx, 64);

  l1c_int i=0;
  // Leave the input vector unchanged.
  double x0_tmp = x_a[0];
  x_a[0] = x_a[0] * sqrt(2.0);

  for (i=0; i<dctmkl_Nx; i++){
    tmp_x[i] = x_a[i];
  }
  tmp_x[dctmkl_Nx] = 0.0; //Strange mkl requirement. dctmkl_tmp_x should have size N+1


  dctmkl_ipar[10] = 1; //Dont normalize
  d_forward_trig_transform(tmp_x, &mkl_dct_handle, dctmkl_ipar, dctmkl_dpar, &dctmkl_stat);

  // Apply the subsampling operation.
  for (i=0; i<dctmkl_Ny; i++){
    y_a[i] = tmp_x[ dctmkl_pix_mask_idx[i]] * dctmkl_root_1_by_2N;
  }

  /*Undo the scaling operation. */
  x_a[0] = x0_tmp;

}


void dctmkl_MtEty( double * restrict y, double * restrict x){
  /*Apply x = M^T * E^T * y.
    -- y should have dimension at least dctmkl_Ny
    -- x should have dimension at least dctmkl_Nx
  */

  l1c_int i=0;
  double *x_a = __builtin_assume_aligned(x, 64);
  double *y_a = __builtin_assume_aligned(y, 64);
  double *tmp_x_a = __builtin_assume_aligned(dctmkl_tmp_x, 64);

  // E^T * y --> y_sparse. y_sparse should be set to all zeros, and pix_mask does not change.
  // Returns a pointer to the array containing the result, which has length N.
  for(i=0; i<dctmkl_Nx+1; i++){
   tmp_x_a[i] = 0.0;
  }

  for (i=0; i<dctmkl_Ny; i++){
    tmp_x_a[dctmkl_pix_mask_idx[i]] = y_a[i];
  }

  // fftw_execute(dctmkl_plan_MtEty); // This is the M^T * dctmkl_y_sparse
  dctmkl_ipar[10] = 1; //Dont normalize
  d_backward_trig_transform(tmp_x_a, &mkl_dct_handle, dctmkl_ipar, dctmkl_dpar, &dctmkl_stat);

  tmp_x_a[0] = tmp_x_a[0]/sqrt(2.0);
  for (i=0; i<dctmkl_Nx; i++){
    x_a[i] = tmp_x_a[i] * dctmkl_root_1_by_2N;;
  }
}


void dctmkl_MtEt_EMx_new(double *x_fftw, double *z){
  /* Performs the Multiplication z = (EM)^T * (EM) * x
     on the supplied vector x.

     -- x should have dimension at least dctmkl_Nx.
     -- z should have dimension at least dctmkl_Nx
  */

  dctmkl_EMx_new(x_fftw, dctmkl_tmp_y); //writes output y to global, which belongs to fftw.
  dctmkl_MtEty(dctmkl_tmp_y, z);


}

#else

#include "l1c_common.h"
void dctmkl_setup(l1c_int Nx, l1c_int Ny, l1c_int *pix_mask_idx){
  Nx +=0;
  Ny +=0;
  pix_mask_idx[0] +=0;
  perror("MKL is not enabled. Recompile with USE_MKL=1");
}
#endif
