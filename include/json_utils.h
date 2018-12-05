#ifndef _JSON_UTILS_
#define  _JSON_UTILS_

long get_file_length(FILE *file_ptr);
int extract_json_double_array(cJSON *data_json, char *name, double **x, int *N);
int load_file_as_text(char *fname, char **file_data);
int load_file_to_json(char *fname, cJSON **data_json);
void print_vec(int N, double *x, char *name);

int extract_json_int(cJSON *data_json, char * name, int *N);
int extract_json_double(cJSON *data_json, char * name, double *a);

int extract_json_int_array(cJSON *data_json, char *name, int **x, int *N);
#endif
