#ifndef _CGSOLVE_
#define _CGSOLVE_
#include "config.h"

#include <stddef.h>
#include "l1c_common.h"
/** Struct containing artifacts of the cgsolve routine. */
typedef struct CgResults_{
  double cgres; /**< Residual */
  l1c_int cgiter;   /**< Number of completed iterations. */

} CgResults;

typedef struct CgParams_{
  l1c_int verbose; /**< If 0, print nothing, if >0, print status every verbose-th iteration. */
  l1c_int max_iter;/**< Maximum number of solver iterations.*/
  double tol;  /**< Solver tolerance.*/
} CgParams;


extern int cgsolve(double *x, double *b, l1c_int N,  double *Dwork,
            void(*AX_func)(l1c_int n, double *x, double *b, void *AX_data),
            void *AX_data, CgResults *cg_result, CgParams cg_params);

int cgsolve_diag_precond(double *x, double *b, double *M_inv_diag, l1c_int N, double *Dwork,
                         void(*AX_func)(l1c_int n, double *x, double *b, void *AX_data), void *AX_data,
                         CgResults *cg_result, CgParams cg_params);

extern void Ax(l1c_int n, double *x, double *b, void *AX_data);

extern void Ax_sym(l1c_int n, double *x, double *b, void *AX_data);

#endif
