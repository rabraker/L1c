#include "config.h"
#include "l1c_timing.h"
#include <math.h>
#include <stdio.h>

#include "TV.h"
#include "bregman.h"
#include "cblas.h"
#include "l1c_logging.h"
#include "l1c_math.h"

/*
  This file contains implementations from "The Split Bregman Method
  for L1-Regularized Problems," by Tom Goldstein and Stanley Osher.

  Currently, only the Anistropic TV denoising problem is implemented.

 */

/* Forward declarations */
static inline void breg_shrink1(l1c_int N, double* x, double* d, double gamma);

static inline void breg_mxpy_z(l1c_int N, double* restrict x, double* restrict y, double* z);

static void
breg_anis_jacobi(int n, int m, double* uk_1, double* uk, double* rhs, double* D, double lambda);

static void breg_anis_guass_seidel(int n, int m, double* u, double* rhs, double mu, double lambda);

static void breg_anis_rhs(l1c_int n,
                          l1c_int m,
                          double* f,
                          double* dx,
                          double* bx,
                          double* dy,
                          double* by,
                          double* rhs,
                          double mu,
                          double lambda,
                          double* dwork1,
                          double* dwork2);

static void hess_inv_diag(l1c_int n, l1c_int m, double mu, double lambda, double* D);

/**
   For unit tests, exports a table of function pointers.
 */

BregFuncs breg_get_functions() {
  BregFuncs bfun = {.breg_shrink1 = breg_shrink1,
                    .breg_mxpy_z = breg_mxpy_z,
                    .breg_anis_guass_seidel = breg_anis_guass_seidel,
                    .breg_anis_rhs = breg_anis_rhs,
                    .hess_inv_diag = hess_inv_diag,
                    .breg_anis_jacobi = breg_anis_jacobi};

  return bfun;
}

void breg_shrink1(l1c_int N, double* restrict x, double* restrict d, double gamma) {
  double* x_ = __builtin_assume_aligned(x, DALIGN);
  double* d_ = __builtin_assume_aligned(d, DALIGN);
  int i = 0;
  double sign_x = 0;

#pragma omp parallel for private(sign_x, i)
  for (i = 0; i < N; i++) {
    sign_x = dsign(x_[i]);
    d_[i] = sign_x * max(fabs(x_[i]) - gamma, 0);
  }
}

/**
   Performs the vector operation
   z = -x + y.
   Assumes x, y, and z are aligned on a DALIGN (default: 64) byte boundary.
*/
void breg_mxpy_z(l1c_int N, double* restrict x, double* restrict y, double* z) {
  double* x_ = __builtin_assume_aligned(x, DALIGN);
  double* y_ = __builtin_assume_aligned(y, DALIGN);
  double* z_ = __builtin_assume_aligned(z, DALIGN);

#pragma omp parallel for
  for (int i = 0; i < N; i++) {
    z_[i] = -x_[i] + y_[i];
  }
}

void breg_hess_eval(l1c_int n,
                    l1c_int m,
                    double* x,
                    double* y,
                    double mu,
                    double lambda,
                    double* dwork1,
                    double* dwork2) {

  l1c_DyTDy(n, m, 1.0, x, dwork1);
  l1c_DxTDx(n, m, 1.0, x, dwork2);
  for (int i = 0; i < n * m; i++) {
    y[i] = mu * x[i] + lambda * dwork1[i] + lambda * dwork2[i];
  }
}

/*
  Not currently using this, but I think it could be useful, especially if I decide to try
  Jacobi iteration instead of guass seidel.
  for n=5, m=4,
  D = [2. 3. 3. 2. 3. 4. 4. 3. 3. 4. 4. 3. 3. 4. 4. 3. 2. 3. 3. 2.]
  row
*/
void hess_inv_diag(l1c_int n, l1c_int m, double mu, double lambda, double* D) {
  l1c_int k = 0, j = 0;

  double mu_2lam = 1.0 / (mu + 2 * lambda);
  double mu_3lam = 1.0 / (mu + 3 * lambda);
  double mu_4lam = 1.0 / (mu + 4 * lambda);

  (void)mu;
  (void)lambda;
  // The first m elements.
  D[0] = mu_2lam;
  for (j = 1; j < m - 1; j++) {
    D[k * m + j] = mu_3lam;
  }
  D[m - 1] = mu_2lam;

  // The middle (m-2)*n elements
  for (k = 1; k < n - 1; k++) {
    D[k * m] = mu_3lam;
    for (j = 1; j < m - 1; j++) {
      D[k * m + j] = mu_4lam;
    }
    D[k * m + j] = mu_3lam;
  }

  // The last m elements.
  D[(n - 1) * m] = mu_2lam;
  for (j = (n - 1) * m + 1; j < m * n - 1; j++) {
    D[j] = mu_3lam;
  }
  D[m * n - 1] = mu_2lam;
}

/* ---------------- JACOBI -----------------------------------------------

DxTDx (4 by 3) (w/o diagonal.):

4 by 3:
[[ 0. -1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
[-1.  0. -1.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
[ 0. -1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
[ 0.  0.  0.  0. -1.  0.  0.  0.  0.  0.  0.  0.]
[ 0.  0.  0. -1.  0. -1.  0.  0.  0.  0.  0.  0.]
[ 0.  0.  0.  0. -1.  0.  0.  0.  0.  0.  0.  0.]
[ 0.  0.  0.  0.  0.  0.  0. -1.  0.  0.  0.  0.]
[ 0.  0.  0.  0.  0.  0. -1.  0. -1.  0.  0.  0.]
[ 0.  0.  0.  0.  0.  0.  0. -1.  0.  0.  0.  0.]
[ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0. -1.  0.]
[ 0.  0.  0.  0.  0.  0.  0.  0.  0. -1.  0. -1.]
[ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0. -1.  0.]]

DyTDy (4
*/

void breg_anis_jacobi(
    int n, int m, double* uk, double* dwork, double* rhs, double* D, double lambda) {

  double ui_min_m = 0, ui_p_m = 0;
  int i = 0, Len = n * m, j = 0;
  /* From l1c_DyTDy, with diagonal removed.*/
  // #pragma omp parallel for private(ui_min_m, ui_p_m, i)
  for (i = 0; i < Len; i++) {
    ui_min_m = i < m ? 0 : uk[i - m];
    ui_p_m = i < Len - m ? uk[i + m] : 0;
    dwork[i] = rhs[i] - lambda * (-ui_min_m - ui_p_m);
  }
  /* From l1c_DxTDx, with diagonal removed.*/
  // #pragma omp parallel for
  for (i = 0; i < n; i++) {
    dwork[i * m] += lambda * uk[i * m + 1];
    for (j = i * m + 1; j < (i + 1) * m - 1; j++) {
      dwork[j] += -lambda * (-uk[j - 1] - uk[j + 1]);
    }
    // We now have col=(row+1)*m
    dwork[j] += -lambda * (-uk[j - 1]);
  }

  l1c_dxmuly_z(Len, D, dwork, uk);
}

/*
  According to https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method,
  Guass-seidel needs to be done sequentially, and u is updated in place.
  Thus, the ordering below is very specific (and assumes that the image data
  in u is stored in row major order): we start the top left corner (i,j=0),
  move across the top row, to the top right pixel (i=0, j=m), the do each row, taking
  special account of the first and last column.

  This is similar to the guass-seidel example on page 357 of "Matrix Computations" by Golub
  and Val Loan. Their example is of the discretized Poisson equation: the difference here is
  in the boundary condition.

  For a 4 by 3 image concatenated in row-major order, the Hessian can be decomposed, for the
  purposes of guass-seidel as

  [ 0., -1.,  0., -1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]
  [-1.,  0., -1.,  0., -1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]
  [ 0., -1.,  0.,  0.,  0., -1.,  0.,  0.,  0.,  0.,  0.,  0.]
  [-1.,  0.,  0.,  0., -1.,  0., -1.,  0.,  0.,  0.,  0.,  0.]
  [ 0., -1.,  0., -1.,  0., -1.,  0., -1.,  0.,  0.,  0.,  0.]
  [ 0.,  0., -1.,  0., -1.,  0.,  0.,  0., -1.,  0.,  0.,  0.]
  [ 0.,  0.,  0., -1.,  0.,  0.,  0., -1.,  0., -1.,  0.,  0.] * lambda + mu*I + lambd*D,
  [ 0.,  0.,  0.,  0., -1.,  0., -1.,  0., -1.,  0., -1.,  0.]
  [ 0.,  0.,  0.,  0.,  0., -1.,  0., -1.,  0.,  0.,  0., -1.]
  [ 0.,  0.,  0.,  0.,  0.,  0., -1.,  0.,  0.,  0., -1.,  0.]
  [ 0.,  0.,  0.,  0.,  0.,  0.,  0., -1.,  0., -1.,  0., -1.]
  [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., -1.,  0., -1.,  0.]

  where D = diag([2., 3., 2., 3., 4., 3., 3., 4., 3., 2., 3., 2.]).

  We want to solve

  H*u = b.

  The guass-seidel update is given by the recursion

           1
   u_i =  -----    [ b_i - \sum_{j~=i} (H_ij * u_j)  ]
          H_(ii)

 */
void breg_anis_guass_seidel(int n, int m, double* u, double* rhs, double mu, double lambda) {

  const double one_by_mu_4lam = 1.0 / (mu + 4 * lambda);
  const double one_by_mu_3lam = 1.0 / (mu + 3 * lambda);
  const double one_by_mu_2lam = 1.0 / (mu + 2 * lambda);
  double tmp;
  int j = 0, i = 0;

  /* Top left left. */
  tmp = rhs[i * m + j] + (u[i * m + j + 1] + u[(i + 1) * m + j]) * lambda;
  u[i * m + j] = tmp * one_by_mu_2lam;

  /* The top row. */
  i = 0;
  for (j = 1; j < m - 1; j++) {
    tmp = rhs[i * m + j] + (u[i * m + j + 1] + u[i * m + j - 1] + u[(i + 1) * m + j]) * lambda;
    u[i * m + j] = tmp * one_by_mu_3lam;
  }
  /* The top, right corner.*/
  j = m - 1;
  i = 0;
  tmp = rhs[i * m + j] + (u[i * m + j - 1] + u[(i + 1) * m + j]) * lambda;
  u[i * m + j] = tmp * one_by_mu_2lam;

  /* The central pixels. The element in the first column
     should be done before the middle of the row, and the last column
     should be done after the middle row.
  */
  for (i = 1; i < n - 1; i++) {
    /* The left pixel in row i.*/
    j = 0;
    tmp = rhs[i * m + j] + (u[i * m + j + 1] + u[(i + 1) * m + j] + u[(i - 1) * m + j]) * lambda;
    u[i * m + j] = tmp * one_by_mu_3lam;
    /* Across the i-th row. */
    for (j = 1; j < m - 1; j++) {
      tmp =
          rhs[i * m + j] +
          (u[i * m + j + 1] + u[i * m + j - 1] + u[(i + 1) * m + j] + u[(i - 1) * m + j]) * lambda;
      u[i * m + j] = tmp * one_by_mu_4lam;
    }
    /* Last column in row i. */
    j = m - 1;
    tmp = rhs[i * m + j] + (u[i * m + j - 1] + u[(i + 1) * m + j] + u[(i - 1) * m + j]) * lambda;
    u[i * m + j] = tmp * one_by_mu_3lam;
  }

  /* The bottom left corner.*/
  j = 0;
  i = n - 1;
  tmp = rhs[i * m + j] + (u[i * m + j + 1] + u[(i - 1) * m + j]) * lambda;
  u[i * m + j] = tmp * one_by_mu_2lam;

  /* The bottom row, excluding corners*/
  i = n - 1;
  for (j = 1; j < m - 1; j++) {
    tmp = rhs[i * m + j] + (u[i * m + j + 1] + u[i * m + j - 1] + u[(i - 1) * m + j]) * lambda;
    u[i * m + j] = tmp * one_by_mu_3lam;
  }

  /* The bottom right corner.*/
  j = m - 1;
  i = n - 1;
  tmp = rhs[i * m + j] + (u[i * m + j - 1] + u[(i - 1) * m + j]) * lambda;
  u[i * m + j] = tmp * one_by_mu_2lam;
}

void breg_anis_rhs(l1c_int n,
                   l1c_int m,
                   double* f,
                   double* dx,
                   double* bx,
                   double* dy,
                   double* by,
                   double* rhs,
                   double mu,
                   double lambda,
                   double* dwork1,
                   double* dwork2) {

  l1c_int N = n * m;
  // RHS = mu*f + lam*Delx^T*(dxk-bxk) + lam*Dely^T(dyk - byk)
  breg_mxpy_z(N, by, dy, dwork1);
  l1c_DyT(n, m, lambda, dwork1, dwork2); // dwork2 = lam*Dely^T(dyk - byk)
  breg_mxpy_z(N, bx, dx, dwork1);
  l1c_DxT(n, m, lambda, dwork1, rhs); // rhs = lam*Delx^T(xyk - bxk)

  cblas_daxpy(N, mu, f, 1, dwork2, 1);    // dwork2 = mu*f + lam*Dely^T(dyk - byk)
  cblas_daxpy(N, 1.0, dwork2, 1, rhs, 1); // RHS
}

/**
 * @defgroup bregman Optimizations which use Bregman splitting.
 * @{*/

/**
 * Given an `n` by `m` image `f`, solves the anistropic TV denoising problem
 * \f[
 * \min_{u} ||\nabla_x u||_1 + ||\nabla_y||_1 + 0.5\mu ||u - f||_2
 * \f]
 * using Bregman Splitting. This algorithm was developed by Tom
 * Goldstein and Stanley Osher in @cite goldstein_splitbregman_2009.
 *
 *
 * @param[in] n Number of rows of the image.
 * @param[in] m Number of rows of the image.
 * @param[out] uk An `n*m` array allocated by malloc_doubl(), which will
 * contain the result of the optimization.
 * @param[in] f An `n*m` array, allocated by malloc_doubl(), which contains the noisy image.
 *
 * @param[in] mu See equation above.
 * @param[in] tol The iterations stop when \f$ \frac{||u_{k-1} - u_k||_2}{||u_k||} < tol\f$
 * @param[in] max_iter Maximum number of iterations.
 * @param[in] max_jac_iter Number of iterations to run the linear solver.
 * The original paper sets this to 1.
 */
int l1c_breg_anistropic_TV(l1c_int n,
                           l1c_int m,
                           double* uk,
                           double* f,
                           double mu,
                           double tol,
                           int max_iter,
                           int max_jac_iter) {

  int iter = 0, N = n * m, status = 0;
  double lambda = 2 * mu, dnrm_err = 0;
  double *uk_1 = NULL, *d_x = NULL, *d_y = NULL, *b_x = NULL, *b_y = NULL;
  double *dwork1 = NULL, *dwork2 = NULL, *rhs = NULL;
  double *Dxu_b = NULL, *Dyu_b = NULL;

  uk_1 = l1c_calloc_double(N);
  d_x = l1c_calloc_double(N);
  d_y = l1c_calloc_double(N);
  b_x = l1c_calloc_double(N);
  b_y = l1c_calloc_double(N);
  dwork1 = l1c_calloc_double(N);
  dwork2 = l1c_calloc_double(N);

  /* Convenience handles. Note below that usage is separated.*/
  Dxu_b = dwork1;
  Dyu_b = dwork2;
  rhs = l1c_calloc_double(N);

  if (!uk_1 || !d_x || !d_y || !b_x || !b_y || !dwork1 || !dwork2 || !rhs) {
    status = L1C_OUT_OF_MEMORY;
    goto exit;
  }

  cblas_dcopy(N, f, 1, uk, 1);

  dnrm_err = INFINITY;

  for (iter = 1; iter <= max_iter; iter++) {

    breg_anis_rhs(n, m, f, d_x, b_x, d_y, b_y, rhs, mu, lambda, dwork1, dwork2);

    for (int k = 1; k <= max_jac_iter; k++) {
      breg_anis_guass_seidel(n, m, uk, rhs, mu, lambda);
    }
    dnrm_err = l1c_dnrm2_rel_err(N, uk, uk_1);

    /* Compute Dyu_b = Del_y*u + b. */
    l1c_Dx(n, m, 1.0, uk, Dxu_b);
    cblas_daxpy(N, 1.0, b_x, 1, Dxu_b, 1);
    /*  Dxu_b = Del_x*u + b. */
    l1c_Dy(n, m, 1.0, uk, Dyu_b);
    cblas_daxpy(N, 1.0, b_y, 1, Dyu_b, 1);

    /* Apply shrink operators. */
    breg_shrink1(N, Dxu_b, d_x, 1.0 / lambda);
    breg_shrink1(N, Dyu_b, d_y, 1.0 / lambda);

    /* Bregman update. */
    breg_mxpy_z(N, d_x, Dxu_b, b_x);
    breg_mxpy_z(N, d_y, Dyu_b, b_y);

    cblas_dcopy(N, uk, 1, uk_1, 1);
    if (dnrm_err < tol) {
      break;
    }
  }

  /*Cleanup before exit. */

exit:
  l1c_free_double(uk_1);
  l1c_free_double(d_x);
  l1c_free_double(d_y);
  l1c_free_double(b_x);
  l1c_free_double(b_y);

  l1c_free_double(dwork1);
  l1c_free_double(dwork2);
  l1c_free_double(rhs);

  return status;
}

/** @}*/
