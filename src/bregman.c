#include <math.h>
#include "mkl.h"
#include "l1qc_common.h"
#include "bregman.h"
#include <stdio.h>
/* Forward declarations */
static void breg_Dx(l1c_int n, l1c_int m, double *X, double *Dx);

static void breg_Dy(l1c_int n, l1c_int m, double *X, double *Dy);

static void breg_shrink1(l1c_int N, double *x, double *d, double gamma);

static void breg_mxpy_z(l1c_int N, double * restrict x, double * restrict y, double *z);

BregFuncs breg_get_functions(){
  BregFuncs bfun = {.breg_Dx=breg_Dx,
                    .breg_Dy=breg_Dy,
                    .breg_shrink1=breg_shrink1,
                    .breg_mxpy_z=breg_mxpy_z};

  return bfun;
}

static void breg_Dx(l1c_int n, l1c_int m, double *X, double *Dx){
  int row=0, col=0;
  for (row=0; row<n; row++){
    for (col=row*m; col<(row+1)*m-1; col++){
      Dx[col] = X[col+1] - X[col];
    }
  }

}

static void breg_Dy(l1c_int n, l1c_int m, double *X, double *Dy){
  int row=0, col=0;
  for (row=0; row<n-1; row++){
    for(col = row * m; col < (row+1) * m; col++){
      Dy[col] = X[col+m] - X[col];
    }
  }
}



void breg_shrink1(l1c_int N, double *x, double *d, double gamma){

  double sign_x=0;
  for (int i=0; i<N; i++){
    sign_x = dsign(x[i]);
    d[i] = sign_x * max(fabs(x[i]) - gamma, 0);
  }
}


/**
   Performs the vector operation
   z = -x + y.
   Assumes x, y, and z are aligned on a DALIGN (default: 64) byte boundary.
*/
void breg_mxpy_z(l1c_int N, double * restrict x, double * restrict y, double *z){
  double *x_ = __builtin_assume_aligned(x, DALIGN);
  double *y_ = __builtin_assume_aligned(y, DALIGN);
  double *z_ = __builtin_assume_aligned(z, DALIGN);

  for (int i=0; i<N; i++){

    z_[i] = -x_[i] + y_[i];
  }

}

/**
   Computes the 2-norm, nrm, of the error between x and y, i.e.,
      nrm = ||x - y||_2
   Assumes x and y are aligned to a DALIGN (default: 64) byte boundary.
 */
double l1c_norm2_err(l1c_int N, double * restrict x, double * restrict y){
  double *x_ = __builtin_assume_aligned(x, DALIGN);
  double *y_ = __builtin_assume_aligned(y, DALIGN);

  double nrm = 0.0, diff = 0.0;
  for (int i=0; i<N; i++){
    diff = x_[i] - y_[i];
    nrm += diff * diff;
  }
  return sqrt(nrm);
}

// void breg_GaussSeidel(l1c_int n, l1c_int m, double *uk_1, double *uk, double *d_x, double *d_y,
//                       double *b_x, double *b_y, double *f, double mu, double lambda){
//   /*
//     In the authors original code, have 5 separate loops to handle the 4 edges and one for the
//     middle of the matrix. Here, I assume that space has been allocated at least one slot before
//     the address of the pointers.
//    */
//   l1c_int row=0, col=0;
//     for (row=0; row<n; row++){
//       for (col=row*m; col<(row+1)*m - 1; col++){
//         Dx[col] = X[col+1] - X[col];
//       }
//     }

//   }

// }


int breg_anistropic_TV(l1c_int n, l1c_int m, double *uk, double lambda, double tol, l1c_int max_iter){
  int iter=0, N=n*m, status=0;

  double *uk_1=NULL, *d_x=NULL, *d_y = NULL, *b_x=NULL, *b_y=NULL;
  double *Dxu_b=NULL, *Dyu_b=NULL;

  uk_1 = malloc_double(N);
  d_x = malloc_double(N);
  d_y = malloc_double(N);
  b_x = malloc_double(N);
  b_y = malloc_double(N);
  Dxu_b = malloc_double(N);
  Dyu_b = malloc_double(N);
  if (!uk_1 | !d_x | !d_y | !b_x | !b_y | !Dxu_b | !Dyu_b){
    status = 1;
  goto exit;
  }

  for (iter=1; iter<=max_iter; iter++){

    //uk_1 = G^k

    /*Compute Dxu_b = Del_x*u + b */
    /*Compute Dyu_b = Del_y*u + b */
    breg_Dx(n, m, uk, Dxu_b);
    breg_Dy(n, m, uk, Dyu_b);
    cblas_daxpy(N, 1.0, b_x, 1, Dxu_b, 1);
    cblas_daxpy(N, 1.0, b_y, 1, Dyu_b, 1);

    /*Apply shrink operators */
    breg_shrink1(N, Dxu_b, d_x, 1.0/lambda);
    breg_shrink1(N, Dyu_b, d_y, 1.0/lambda);

    /*Bregman update */
    breg_mxpy_z(N, d_x, Dxu_b, b_x);
    breg_mxpy_z(N, d_y, Dyu_b, b_y);


    if (l1c_norm2_err(N, uk, uk_1) < tol){
      cblas_dcopy(N, uk, 1, uk_1, 1);
      break;
    }
    cblas_dcopy(N, uk, 1, uk_1, 1);
  }

  /*Cleanup before exit. */
  goto exit;

 exit:
  free_double(uk_1);
  free_double(d_x);
  free_double(d_y);
  free_double(b_x);
  free_double(b_y);
  free_double(Dxu_b);
  free_double(Dyu_b);

  return status;
}
