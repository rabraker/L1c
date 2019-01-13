
#ifndef _CHECK_UTILS_
#define _CHECK_UTILS_
#include <check.h>
#include <math.h>


/* Floating point tolerance comparison macros with improved output
 * compared to ck_assert(). */
/* OP can have values: <; >=. */
#define _ck_assert_floating_array_absdiff_op_tol(N, X, Y, OP, T, TP, TM) \
  do {                                                                  \
  for (int i=0; i<N;i++){                                               \
    TP _ck_x = (X[i]);                                                  \
    TP _ck_y = (Y[i]);                                                  \
    TP _ck_t = (T);                                                     \
    ck_assert_msg(fabsl(_ck_y - _ck_x) OP _ck_t,                        \
                  "\n     Assertion '%s' failed: %s[%d] == %.*" TM "g, %s[%d] == %.*" TM "g, %s == %.*" TM "g", \
                  "fabsl("#Y" - "#X") "#OP" "#T,                        \
                  #X, i, (int)CK_FLOATING_DIG, _ck_x,                    \
                  #Y, i, (int)CK_FLOATING_DIG, _ck_y,                    \
                  #T, (int)CK_FLOATING_DIG, _ck_t);                     \
    }                                                                   \
  } while (0)

/* Floating point tolerance comparison macros with improved output
 * compared to ck_assert(). */
/* OP can have values: <; >=. */
#define _ck_assert_floating_array_absdiff_op_reltol(N, X, Y, OP, T, TP, TM) \
  do {                                                                  \
    TP _ck_nrm_x=0, _ck_nrm_y=0;                                        \
    for (int i=0; i<N; i++){                                           \
      _ck_nrm_x += X[i] * X[i];                                         \
      _ck_nrm_y += Y[i] * Y[i];                                         \
    }                                                                   \
    _ck_nrm_x = (TP) sqrt((double)_ck_nrm_x);                           \
    _ck_nrm_y = (TP) sqrt((double)_ck_nrm_y);                           \
                                                                        \
  for (int i=0; i<N;i++){                                               \
    TP _ck_x = (X[i]);                                                  \
    TP _ck_y = (Y[i]);                                                  \
    TP _ck_t = (T);                                                     \
    ck_assert_msg(fabsl(_ck_y/_ck_nrm_y - _ck_x/_ck_nrm_x) OP _ck_t,    \
                  "\n     Assertion '%s' failed: %s[%d] == %.*" TM "g, %s[%d] == %.*" TM "g, %s == %.*" TM "g", \
                  "fabsl("#Y" - "#X") "#OP" "#T,                        \
                  #X, i, (int)CK_FLOATING_DIG, _ck_x,                    \
                  #Y, i, (int)CK_FLOATING_DIG, _ck_y,                    \
                  #T, (int)CK_FLOATING_DIG, _ck_t);                     \
    }                                                                   \
  } while (0)


#define ck_assert_double_array_eq_reltol(N, X, Y, T)  _ck_assert_floating_array_absdiff_op_reltol(N, X, Y, <, T, double, "")
/**
 * Check two double precision floating point array to determine if not X[i]/nrm≈Y[i]/nrm
 * with specified tolerance, for i=0,...N-1
 *
 * If X ≈ Y with error < T, the test fails.
 *
 * @param X floating point number (double)
 * @param Y floating point number (double) to compare against X
 * @param T tolerance (double)
 *
 * @note If the check fails, the remaining of the test is aborted
 *
 * @since 0.11.0
 */



#define ck_assert_double_array_eq_tol(N, X, Y, T)  _ck_assert_floating_array_absdiff_op_tol(N, X, Y, <, T, double, "")
/**
 * Check two double precision floating point array to determine if not X[i]/norm≈Y[i]
 * with specified tolerance, for i=0,...N-1
 *
 * If X ≈ Y with error < T, the test fails.
 *
 * @param X floating point number (double)
 * @param Y floating point number (double) to compare against X
 * @param T tolerance (double)
 *
 * @note If the check fails, the remaining of the test is aborted
 *
 * @since 0.11.0
 */


#endif
