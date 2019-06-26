#include <cblas.h>
#include <math.h>

#include "l1c.h"
#include "l1qc_newton.h"
#include "l1c_math.h"
#include "linesearch.h"
#include "l1c_logging.h"

/**
   Performs a backtracking linesearch for objective functions with (x, u) variables.

   Searches for the largest `s` such that
   \f[
   f_1(x + s*dx) + f_2(u + s*du)< f_1(x) + f_2(u) + alpha * [\nabla_x, \nabla_u]^T [dx\\du]
   \f]

   The algorithm starts at `s = ls_params.s` and successively decreases `s` by the
   factor `ls_params.beta`. The algorithm will terminate once either the line search succeeds,
   of the number of iterations surpasses `ls_params.max_iter`. If the number of iterations
   surpasses `ls_params.max_iter`, the line search has failed and the routine returns
   with `ls_stat.status = 1`.


   @param[in] m Length of the vectors `x`.

   @param[in,out] x On entry, contains the current value of the optimization variable. If the
   linesearch was succesfull, then on exit `x = x + s*dx`, where is `s` is the step size which
   caused the line search to succeed. Assumed to by aligned to DALIGN.

   @param[in,out] u. On entry, contains the current value of the optimization variable. If the
   linesearch was succesfull, then on exit `u = u + s*du`, where is `s` is the step size which
   caused the line search to succeed. Assumed to by aligned to DALIGN.

   @param[in] dx The x part of decent direction. Assumed to by aligned to DALIGN.

   @param[in] du The u part of decent direction. Assumed to by aligned to DALIGN.

   @param[in,out] fx0 (scalar) holds the current value of the objective function.
   If the line search succeeds, then on exit, will contain the objective function
   value evalute at x + s*dx.

   @param[in,out] g_dot_dxu (scalar) value [\nabla_x, \nabla_y] [dx; du]

   @param[in,out] prob_data Any auxilary data needed by `objective_fun`. Some objective functions
   will leave computed auxilary variables in obj_data. If this is the case, then on exit,
   obj_data will contain data associated with the last call to objective_fun. This will be
   the case wheather or not the line search succeeds, because we do not know the structure
   of the data, and thus dont know how to copy it. It is therefore the responsibilty of the
   caller to handle this. In particular, if the line search fails, then `obj_data` may well
   contain garbage (e.g., NAN or INF).

   @param[in] feval A pointer to a function which will compute the value
   of the objective function.

   @param[in] ls_params Struct of line search parameters.

   @param[in] DWORK2 2D Work array with two arrays of size m.

 */

LSStat l1c_linesearch_xu(l1c_int m, double *x, double *u, double *dx, double *du,
                             double *fx0, double g_dot_dxu, void *prob_data,
                             double (*feval)(void *data, double *x, double *u),
                             LSParams ls_params, double **DWORK2){
  LSStat ls_stat = {.flin=0, .step=0, .status=0};
  int iter=0;
  double flin = 0.0, f_val=0.0;
  double fx0_orig=*fx0;

  double step = ls_params.s;
  /* Make things accessible */
  double *xp = DWORK2[0];
  double *up = DWORK2[1];

  for (iter=1; iter<=MAX_LINESEARCH_ITER; iter++){
    /* xp = x + s*dx etc*/
    l1c_daxpy_z(m, step, dx, x, xp);
    l1c_daxpy_z(m, step, du, u, up);

    /* Evaluate functional at the new point. */
    f_val = feval(prob_data, xp, up);

    if ( isnan(f_val)) {
      step = ls_params.beta * step;
      continue;
    }

    /*
      De-rated linear approximation at the current point. We have
      g_dot_dxu = gradf'*[dx; du]
      So flin = f(x0) + alp * s * [Del x; Del u] * [dx, du]
    */
    flin = fx0_orig + ls_params.alpha * step * g_dot_dxu;

    if (f_val <= flin){ /* Sufficient decrease test */
      cblas_dcopy(m, xp, 1, x, 1);
      cblas_dcopy(m, up, 1, u, 1);
      *fx0 = f_val;
      ls_stat.flin = flin;
      ls_stat.step = step;
      ls_stat.iter = iter;
      ls_stat.status = 0;
      return ls_stat;
    }

    step = ls_params.beta * step;
  }

  char spc[] = "   ";
  l1c_printf("%sBacktracking line search failed, returning previous iterate.\n", spc);
  l1c_printf("%sLast line-search values:\n", spc);
  l1c_printf("%s                        iter = %d\n", spc, iter);
  l1c_printf("%s                        fp   = %.10e\n", spc, f_val);
  l1c_printf("%s                        flin = %.10e\n", spc, flin);
  l1c_printf("%s                        fp-flin = %.10f\n", spc, f_val - flin);
  l1c_printf("%s                        step = %.10f\n", spc, step);
  l1c_printf("%s                        ls_params.s = %.10f\n", spc, ls_params.s);

  ls_stat.flin = flin;
  ls_stat.step = step;
  ls_stat.iter = min(iter, MAX_LINESEARCH_ITER);
  ls_stat.status = 1;
  return ls_stat;
}


/**
   Performs a backtracking linesearch.

   Searches for the largest `s` such that
   \f[
   f(x + s*dx) < f(x) + alpha * (\nabla f)^T dx
   \f]

   The algorithm starts at `s = ls_params.s` and successively decreases `s` by the
   factor `ls_params.beta`. The algorithm will terminate once either the line search succeeds,
   of the number of iterations surpasses `ls_params.max_iter`. If the number of iterations
   surpasses `ls_params.max_iter`, the line search has failed and the routine returns
   with `ls_stat.status = 1`.

   @param[in,out] fx0 (scalar) holds the current value of the objective function.
   If the line search succeeds, then on exit, will contain the objective function
   value evalute at x + s*dx.

   @param[in] Length of the vectors `x`.

   @param[in,out]. On entry, contains the current value of the optimization variable. If the
   linesearch was succesfull, then on exit `x = x + s*dx`, where is `s` is the step size which
   caused the line search to succeed. Assumed to by aligned to DALIGN.

   @param[in] dx The decent direction. Assumed to by aligned to DALIGN.

   @param[in] ls_params Struct of line search parameters.

   @param[in] DWORK A work array of length `n`.  Assumed to by aligned to DALIGN.

   @param[in] objective_fun A pointer to a function which will compute the value
   of the objective function.

   @param[in,out] obj_data Any auxilary data needed by `objective_fun`. Some objective functions
   will leave computed auxilary variables in obj_data. If this is the case, then on exit,
   obj_data will contain data associated with the last call to objective_fun. This will be
   the case wheather or not the line search succeeds, because we do not know the structure
   of the data, and thus dont know how to copy it. It is therefore the responsibilty of the
   caller to handle this. In particular, if the line search fails, then `obj_data` may well
   contain garbage (e.g., NAN or INF).
 */
LSStat
l1c_linesearch(l1c_int n, double *x, double *dx, double *fx0, double g_dot_dx, void *obj_data,
               double (*objective_fun)(void *data, double *x), LSParams ls_params, double *DWORK){


  LSStat ls_stat = {.flin=0, .step=0, .status=0};
  int iter = 0;
  double flin = 0.0, fx0_prime = 0.0;

  double step = ls_params.s;

  /* Make things accessible */
  double *xp = DWORK;

  for (iter=1; iter<=MAX_LINESEARCH_ITER; iter++){
    /* xp = x + s*dx */
    l1c_daxpy_z(n, step, dx, x, xp);

    /*Re-evaluate functional. */
    fx0_prime = objective_fun(obj_data, xp);

    if ( isnan(fx0_prime)) {
      step = ls_params.beta * step;
      continue;
    }

    /* need gradf'*dx
     */
    flin = *fx0 + ls_params.alpha * step * g_dot_dx;

    if (fx0_prime <= flin){ /* Sufficient decrease test */
      cblas_dcopy(n, xp, 1, x, 1);
      *fx0 = fx0_prime;

      ls_stat.flin = flin;
      ls_stat.step = step;
      ls_stat.iter = iter;
      ls_stat.status = 0;
      return ls_stat;
    }

    step = ls_params.beta * step;
  }

  char spc[] = "   ";
  l1c_printf("%sBacktracking line search failed, returning previous iterate.\n", spc);
  l1c_printf("%sLast line-search values:\n", spc);
  l1c_printf("%s                        iter = %d\n", spc, iter);
  l1c_printf("%s                        fp   = %.10e\n", spc, fx0_prime);
  l1c_printf("%s                        flin = %.10e\n", spc, flin);
  l1c_printf("%s                        fp-flin = %.10f\n", spc, fx0_prime - flin);
  l1c_printf("%s                        step = %.10f\n", spc, step);
  l1c_printf("%s                        ls_params.s = %.10f\n", spc, ls_params.s);

  ls_stat.flin = flin;
  ls_stat.step = step;
  ls_stat.iter = min(iter, MAX_LINESEARCH_ITER);
  ls_stat.status = 1;
  return ls_stat;
}
