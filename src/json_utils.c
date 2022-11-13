#include "config.h"

#include <cjson/cJSON.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "json_utils.h"
#include "l1c.h"

/* Utility functions */

void print_vec(l1c_int N, double* x, char* name) {
  l1c_int i = 0;
  for (i = 0; i < N; i++) {
    printf("%s[%i] = %f\n", name, (int)i, x[i]);
  }
}

int load_file_as_text(char* fname, char** file_data) {
  int file_len;
  FILE* file_ptr;

  file_ptr = fopen(fname, "r"); /* Open in TEXT mode */

  if (file_ptr == NULL) { /* Could not open file */
    fprintf(stderr, "Error opening %s: %s (%u)\n", fname, strerror(errno), errno);
    return 1;
  }

  file_len = get_file_length(file_ptr);
  *file_data = calloc(file_len + 1, sizeof(char));

  if (*file_data == NULL) {
    fprintf(stderr, "\nInsufficient memory to read file.\n");
    fclose(file_ptr);
    return 1;
  }

  /* Read the entire file into cFile */
  const size_t nmemb = 1;
  const size_t num_read = fread(*file_data, file_len, nmemb, file_ptr);
  fclose(file_ptr);
  if (nmemb == num_read) {
    return 0;
  } else {
    return 1;
  }
}

int load_file_to_json(char* fname, cJSON** data_json) {
  char* file_data = NULL;

  if (load_file_as_text(fname, &file_data))
    return 1;

  *data_json = cJSON_Parse(file_data);

  free(file_data);

  if (data_json == NULL) {
    const char* error_ptr = cJSON_GetErrorPtr();
    if (error_ptr != NULL) {
      fprintf(stderr, "Error before: %s\n", error_ptr);
      return 1;
    }
  }

  return 0;
}

int extract_json_int(cJSON* data_json, char* name, l1c_int* N) {
  if (!data_json) {
    fprintf(stderr, "Error in extract_json_int: cannot extract %s\n", name);
    fprintf(stderr, "         cJSON input is NULL\n");
    return 1;
  }

  cJSON* N_json = cJSON_GetObjectItemCaseSensitive(data_json, name);
  if (!cJSON_IsNumber(N_json)) {
    fprintf(stderr, "Error in extract_json_int: cannot extract %s\n", name);
    return 1;
  }

  *N = (*N_json).valueint;
  return 0;
}

int extract_json_double(cJSON* data_json, char* name, double* a) {
  if (!data_json) {
    fprintf(stderr, "Error in extract_json_double: cannot extract %s\n", name);
    fprintf(stderr, "         cJSON input is NULL\n");
    return 1;
  }

  cJSON* a_json = cJSON_GetObjectItemCaseSensitive(data_json, name);
  if (!cJSON_IsNumber(a_json)) {
    fprintf(stderr, "Error in extract_json_double: cannot extract %s\n", name);
    return 1;
  }

  *a = (*a_json).valuedouble;
  return 0;
}

int extract_json_double_array(cJSON* data_json, char* name, double** x, l1c_int* N) {
  if (!data_json) {
    fprintf(stderr, "Error in extract_json_double_array: cannot extract %s\n", name);
    fprintf(stderr, "         cJSON input is NULL\n");
    return 1;
  }

  int status = 0;
  cJSON* x_json = cJSON_GetObjectItemCaseSensitive(data_json, name);
  if (!cJSON_IsArray(x_json)) {
    fprintf(stderr, "Error in extract_json_double_array: cannot extract %s\n", name);
    status = 1;
    *N = 0;
    goto end;
  }
  *N = cJSON_GetArraySize(x_json);
  *x = l1c_malloc_double(*N);

  if (!*x) {
    status = 1;
    fprintf(stderr, "error allocating memory\n");
    goto end;
  }

  int k = 0;
  cJSON* xk_json;
  cJSON_ArrayForEach(xk_json, x_json) {
    // note to self: if you do *x[k], that would be
    // the next pointer in an array of pointers.
    (*x)[k] = xk_json->valuedouble;
    k = k + 1;
  }

end:
  return status;
}

int extract_json_int_array(cJSON* data_json, char* name, l1c_int** x, l1c_int* N) {
  if (!data_json) {
    fprintf(stderr, "Error in extract_json_int_array: cannot extract %s\n", name);
    fprintf(stderr, "         cJSON input is NULL\n");
    return 1;
  }

  int status = 0;
  cJSON* x_json = cJSON_GetObjectItemCaseSensitive(data_json, name);
  if (!cJSON_IsArray(x_json)) {
    fprintf(stderr, "Error in extract_json_int_array: cannot extract %s\n", name);
    status = 1;
    *N = 0;
    goto end;
  }
  *N = cJSON_GetArraySize(x_json);
  *x = calloc(*N, sizeof(l1c_int));

  int k = 0;
  cJSON* xk_json;
  cJSON_ArrayForEach(xk_json, x_json) {
    // note to self: if you do *x[k], that would be
    // the next pointer in an array of pointers.
    (*x)[k] = xk_json->valueint;
    k = k + 1;
  }

end:
  return status;
}

long get_file_length(FILE* file_ptr) {
  // Will return file at begining of file.
  fseek(file_ptr, 0L, SEEK_END);   /* Position to end of file */
  long file_len = ftell(file_ptr); /* Get file length */
  rewind(file_ptr);                /* Back to start of file */

  return file_len;
}
