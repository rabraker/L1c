#include <cblas.h>
#include <stdio.h>

#include "cblas.h"

/* This function will execute when libl1c.so is loaded, and warn the user if
   it detects that openblas is configured to use pthread instead of openmp
 */
__attribute__((constructor)) void check_opennblas_pthreads() {
  if (openblas_get_parallel() == OPENBLAS_THREAD) {
    fprintf(stderr,
            "***************************************************************\n");
    fprintf(stderr, "                    WARNING \n");
    fprintf(stderr,
            "Openblas is configured with pthreads, not openmp.\n "
            "This is likely to be extremely slow, and interfere with fftw\n"
            "You can inst all the openmp version with \n"
            "apt instal libopenblas-openmp-dev \n"
            "or mitigate the issue (with some performance drop) with\n"
            "OMP_NUM_THREADS=2\n"
            "See also https://github.com/xianyi/OpenBLAS/issues/3187 \n"
            "\n\n");
    fprintf(stderr, "**************************************************************\n");
  }
}
