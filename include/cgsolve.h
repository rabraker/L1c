#ifndef _CGSOLVE_
#define _CGSOLVE_

#include <stddef.h>
/** Struct containing artifacts of the cgsolve routine. */
typedef struct CgResults_{
  double cgres; /**< Residual */
  int cgiter;   /**< Number of completed iterations. */

} CgResults;

typedef struct CgParams_{
  int verbose; /**< If 0, print nothing, if >0, print status every verbose-th iteration. */
  int max_iter;/**< Maximum number of solver iterations.*/
  double tol;  /**< Solver tolerance.*/
} CgParams;


extern int cgsolve(double *x, double *b, size_t n_b,  double *Dwork,
            void(*AX_func)(int n, double *x, double *b, void *AX_data),
            void *AX_data, CgResults *cg_result, CgParams cg_params);


extern void Ax(int n, double *x, double *b, void *AX_data);

extern void Ax_sym(int n, double *x, double *b, void *AX_data);

#endif
