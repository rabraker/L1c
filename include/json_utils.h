#ifndef _JSON_UTILS_
#define  _JSON_UTILS_
#include "l1qc_common.h"

void free_json_text_data(void);

long get_file_length(FILE *file_ptr);

int extract_json_double_array(cJSON *data_json, char *name, double **x, l1c_int *N);

int extract_json_double_array_fftw(cJSON *data_json, char *name, double **x, l1c_int *N);

int load_file_as_text(char *fname, char **file_data);

int load_file_to_json(char *fname, cJSON **data_json);

void print_vec(l1c_int N, double *x, char *name);

int extract_json_int(cJSON *data_json, char * name, l1c_int *N);

int extract_json_double(cJSON *data_json, char * name, double *a);

int extract_json_int_array(cJSON *data_json, char *name, l1c_int **x, l1c_int *N);


#endif
