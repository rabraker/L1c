#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cjson/cJSON.h>

#include "json_utils.h"
#include "l1c.h"

int main(int argc, char** argv) {

  int verbose = 0;
  char* fpath;
  char fpath_[] = "example_img_data.json";

  if (argc >= 2) {
    fpath = malloc(strlen(argv[1]) * sizeof(char) + 1);
    if (!fpath) {
      fprintf(stderr, "Error allocating memory\n");
      return 1;
    }
    strcpy(fpath, argv[1]);
  } else {
    fpath = malloc(strlen(fpath_) * sizeof(char) + 1);
    if (!fpath) {
      fprintf(stderr, "Error allocating memory\n");
      return 1;
    }
    strcpy(fpath, fpath_);
  }

  if (argc >= 3) {
    verbose = atoi(argv[2]);
  } else {
    verbose = 0;
  }

  cJSON* test_data_json = NULL;
  double *x = NULL, *b = NULL;
  l1c_int* pix_idx = NULL;
  l1c_int m = 0, n = 0, status = 0;

  l1c_LBResult lb_res;
  l1c_L1qcOpts l1qc_opts = {.epsilon = .1,
                            .mu = 10,
                            .lbtol = 1e-3,
                            .newton_tol = 1e-3,
                            .newton_max_iter = 50,
                            .verbose = verbose,
                            .l1_tol = 1e-5,
                            .cg_tol = 1e-8,
                            .cg_maxiter = 200,
                            .cg_verbose = 0,
                            .warm_start_cg = 0,
                            .dct_mode = dct1};

  if (load_file_to_json(fpath, &test_data_json)) {
    fprintf(stderr, "Error loading data in l1qc_dct_c\n");
    status = L1C_FILE_READ_FAILURE;
    goto exit;
  }

  status += extract_json_double_array(test_data_json, "b", &b, &n);
  status += extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &n);
  status += extract_json_int(test_data_json, "mtot", &m);

  x = l1c_malloc_double(m);

  if (status || !x) {
    fprintf(stderr, "Error reading JSON file or allocating memory\n");
    status += 1;
    goto exit;
  }
  // time_t start = time(NULL);

  l1qc_dct(m, 1, x, n, b, pix_idx, l1qc_opts, &lb_res);

exit:
  free(fpath);
  l1c_free_double(x);
  free(pix_idx);
  l1c_free_double(b);
  cJSON_Delete(test_data_json);

  return status;
}
