#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "l1c_common.h"
#include "l1qc_newton.h"
#include "omp.h"
#include <cjson/cJSON.h>
#include "json_utils.h"

/* dct_mkl.h defines the Ax and Aty operations.
   To adapt this mex file to a different set of transformations,
   this largely what must be changed. Your new set of transformations
   must expose three functions with the following prototype:

   void Ax(double *x, double *y)
   void Aty(double *y, double *x)
   void AtAx(double *x, double *z)

   where x and z have length n, and y has length m and m<n.
*/

#include "dct_mkl.h"


typedef struct L1qcDctOpts{
  double epsilon;
  double mu;
  double lbtol;
  double tau;
  double lbiter;
  double newton_tol;
  int newton_max_iter;
  int verbose;
  double l1_tol;
  double cgtol;
  int cgmaxiter;
  int warm_start_cg;
}L1qcDctOpts;

int l1qc_dct(int N, double *x_out, int M, double *b, l1c_int *pix_idx,
             L1qcDctOpts opts, LBResult *lb_res);



int main(){

  mkl_set_interface_layer(MKL_INTERFACE_ILP64);
  mkl_set_threading_layer(MKL_THREADING_GNU);

  cJSON *test_data_json;
  double *x=NULL, *b=NULL;
  l1c_int *pix_idx=NULL;
  l1c_int N=0, M=0, status=0;
  char fpath[] = "example_img_data.json";

  LBResult lb_res;
  L1qcDctOpts l1qc_dct_opts = {.epsilon=.1,
							   .mu = 10,
							   .lbtol = 1e-3,
							   .newton_tol = 1e-3,
                               .newton_max_iter = 50,
                               .verbose = 2,
                               .l1_tol = 1e-5,
                               .cgtol = 1e-8,
                               .cgmaxiter = 200,
                               .warm_start_cg=0};


  if (load_file_to_json(fpath, &test_data_json)){
    fprintf(stderr, "Error loading data in test_get_gradient\n");
    goto exit;
  }

  status +=extract_json_double_array(test_data_json, "b", &b, &M);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &M);
  status +=extract_json_int(test_data_json, "N", &N);

  x = malloc_double(N);

  if (status || !x){
    fprintf(stderr, "Error reading JSON file or allocating memory\n");
    goto exit;
  }
  // time_t start = time(NULL);

  l1qc_dct(N, x, M, b, pix_idx, l1qc_dct_opts, &lb_res);

 exit:
  free_double(x);
  free(pix_idx);
  free_double(b);

}
