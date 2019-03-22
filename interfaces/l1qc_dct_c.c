#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "omp.h"
#include <cjson/cJSON.h>

#include "json_utils.h"
#include "l1c_memory.h"
#include "l1c.h"


int main(int argc, char **argv){
#if defined(HAVE_LIBMKL_RT)
  mkl_set_interface_layer(MKL_INTERFACE_ILP64);
  mkl_set_threading_layer(MKL_THREADING_GNU);
#endif
  char *fpath;
  char fpath_[] = "example_img_data.json";

  if (argc >= 2){
    fpath = malloc(strlen(argv[1])*sizeof(char)+1);
    if (!fpath){
      fprintf(stderr, "Error allocating memory\n");
      return 1;
    }
    strcpy(fpath, argv[1]);
  }else{
    fpath = malloc(strlen(fpath_)*sizeof(char)+1);
    if (!fpath){
      fprintf(stderr, "Error allocating memory\n");
      return 1;
    }
    strcpy(fpath, fpath_);
  }

  cJSON *test_data_json=NULL;
  double *x=NULL, *b=NULL;
  l1c_int *pix_idx=NULL;
  l1c_int N=0, M=0, status=0;



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
    fprintf(stderr, "Error loading data in l1qc_dct_c\n");
    status = L1C_FILE_READ_FAILURE;
    goto exit;
  }

  status +=extract_json_double_array(test_data_json, "b", &b, &M);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &M);
  status +=extract_json_int(test_data_json, "N", &N);

  x = malloc_double(N);

  if (status || !x){
    fprintf(stderr, "Error reading JSON file or allocating memory\n");
    status +=1;
    goto exit;
  }
  // time_t start = time(NULL);

  l1qc_dct(N, 1, x, M, b, pix_idx, l1qc_dct_opts, &lb_res);

 exit:
  free(fpath);
  free_double(x);
  free(pix_idx);
  free_double(b);
  cJSON_Delete(test_data_json);

  return status;

}
