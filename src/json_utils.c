#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <cjson/cJSON.h>

#include "json_utils.h"
#include "l1c_common.h"


/* Utility functions */


void print_vec(l1c_int N, double *x, char *name){
  l1c_int i = 0;
  for (i=0; i<N; i++){
    printf("%s[%i] = %f\n", name, (int)i, x[i]);
  }
}


int load_file_as_text(char *fname, char **file_data){
  int file_len;
  FILE *file_ptr;

  file_ptr = fopen(fname, "r");  /* Open in TEXT mode */

  if (file_ptr == NULL ){             /* Could not open file */
    printf("Error opening %s: %s (%u)\n", fname, strerror(errno), errno);
    return 1;
  }

  file_len = get_file_length(file_ptr);
  *file_data = calloc(file_len + 1, sizeof(char));

  if(*file_data == NULL ){
    printf("\nInsufficient memory to read file.\n");
    fclose(file_ptr);
    return 1;
  }

  fread(*file_data, file_len, 1, file_ptr); /* Read the entire file into cFile */
  fclose(file_ptr);                    /* Close it */

  return 0;
}

int load_file_to_json(char *fname, cJSON **data_json){
  char *file_data =NULL;

  if (load_file_as_text(fname, &file_data))
    return 1;

  *data_json = cJSON_Parse(file_data);

  free(file_data);

  if (data_json == NULL){
    const char *error_ptr = cJSON_GetErrorPtr();
    if (error_ptr != NULL){
      fprintf(stderr, "Error before: %s\n", error_ptr);
      return 1;
    }
  }

  return 0;
}

int extract_json_int(cJSON *data_json, char * name, l1c_int *N){
  cJSON *N_json = cJSON_GetObjectItemCaseSensitive(data_json, name);
  if (!cJSON_IsNumber(N_json) ){
    return 1;
  }

  *N = (*N_json).valueint;
  return 0;
}

int extract_json_double(cJSON *data_json, char * name, double *a){
  cJSON *a_json = cJSON_GetObjectItemCaseSensitive(data_json, name);
  if (!cJSON_IsNumber(a_json)){
    return 1;
  }

  *a = (*a_json).valuedouble;
  return 0;
}


int extract_json_double_array(cJSON *data_json, char *name, double **x, l1c_int *N){
  int status = 0;
  cJSON *x_json = cJSON_GetObjectItemCaseSensitive(data_json, name);
  if (!cJSON_IsArray(x_json) ){
      status = 1;
      *N = 0;
      goto end;
    }
  *N = cJSON_GetArraySize(x_json);
  // *x = calloc(*N, sizeof(double));
  *x = malloc_double(*N);

  if (!*x){
    status =1;
    perror("error allocating memory\n");
    goto end;
  }

  int k=0;
  cJSON *xk_json;
  cJSON_ArrayForEach(xk_json, x_json){
    // note to self: if you do *x[k], that would be
    // the next pointer in an array of pointers.
    (*x)[k] = xk_json->valuedouble;
    k = k+1;
  }

 end:
  return status;
}

int extract_json_int_array(cJSON *data_json, char *name, l1c_int **x, l1c_int *N){
  int status = 0;
  cJSON *x_json = cJSON_GetObjectItemCaseSensitive(data_json, name);
  if (!cJSON_IsArray(x_json) ){
    status = 1;
    *N = 0;
    goto end;
  }
  *N = cJSON_GetArraySize(x_json);
  *x = calloc(*N, sizeof(l1c_int));

  int k=0;
  cJSON *xk_json;
  cJSON_ArrayForEach(xk_json, x_json){
    // note to self: if you do *x[k], that would be
    // the next pointer in an array of pointers.
    (*x)[k] = xk_json->valueint;
    k = k+1;
  }

 end:
  return status;
}


long get_file_length(FILE *file_ptr){
  // Will return file at begining of file.
  fseek(file_ptr, 0L, SEEK_END);     /* Position to end of file */
  long file_len = ftell(file_ptr);        /* Get file length */
  rewind(file_ptr);                  /* Back to start of file */

  return file_len;
}
