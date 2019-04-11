
#ifndef _CHECK_UTILS_
#define _CHECK_UTILS_
#include <check.h>
#include <math.h>

/* This is only defined in new versions of check.
   Since we provide backward compatibility, we must check
   and define if check.h didnt do it.
 */
#ifndef CK_FLOATING_DIG
#define CK_FLOATING_DIG 10
#endif


static inline double ck_min(double a,double b){
  return a < b ? a : b;
}


#ifndef ck_assert_double_eq_tol
/*
  Simplified version of macro availible in check version > 0.11
 */
#define ck_assert_double_eq_tol(X, Y, T)                                \
  do {                                                                  \
    double _ck_x = (X);                                                 \
    double _ck_y = (Y);                                                 \
    double _ck_t = (T);                                                 \
    ck_assert_msg(fabsl(_ck_y - _ck_x) < _ck_t,                         \
                  "Assertion '%s' failed: %s == %.* g, %s == %.* g, %s == %.* g", \
                  "fabsl("#Y" - "#X") < "#T,                            \
                  #X, (int)CK_FLOATING_DIG, _ck_x,                      \
                  #Y, (int)CK_FLOATING_DIG, _ck_y,                      \
                  #T, (int)CK_FLOATING_DIG, _ck_t);                     \
  } while (0)

#endif

#ifndef _ck_assert_floating
/* Floating point number comparison macros with improved output
 * compared to ck_assert(). */
/* OP may be any comparison operator, TP is type, TM is type modifier. */
#define _ck_assert_double(X, OP, Y) do {                            \
    double _ck_x = (X);                                                 \
    double _ck_y = (Y);                                                 \
    ck_assert_msg(_ck_x OP _ck_y,                                       \
                  "Assertion '%s' failed: %s == %.* g, %s == %.* g", \
                  #X" "#OP" "#Y,                                        \
                  #X, (int)CK_FLOATING_DIG, _ck_x,                      \
                  #Y, (int)CK_FLOATING_DIG, _ck_y);                     \
  } while (0)
#endif

#ifndef ck_assert_double_le
/*
 * Check two double precision floating point numbers to determine if X <= Y
 *
 * Availible in check 0.11.0
 */
#define ck_assert_double_le(X, Y) _ck_assert_double(X, <=, Y)

#endif

/* Floating point tolerance Floating point relative tolerance comparison
   for arrays. Based on ck_assert_floating_absdiff_tol.
   OP can have values: <; >=. */
#define _ck_assert_floating_array_absdiff_op_tol(N, X, Y, OP, T, TP, TM) \
  do {                                                                   \
    for (int _ck_i_=0; _ck_i_<N;_ck_i_++){                               \
    TP _ck_x = (X[_ck_i_]);                                              \
    TP _ck_y = (Y[_ck_i_]);                                              \
    TP _ck_t = (T);                                                      \
    ck_assert_msg(fabsl(_ck_y - _ck_x) OP _ck_t,                         \
                  "\n     Assertion '%s' failed: \n"                     \
                  "      %s[%d] == %.*" TM "g,"                          \
                  " %s[%d] == %.*" TM "g, %s == %.*" TM "g",             \
                  "fabsl("#Y" - "#X") "#OP" "#T,                         \
                  #X, _ck_i_, (int)CK_FLOATING_DIG, _ck_x,               \
                  #Y, _ck_i_, (int)CK_FLOATING_DIG, _ck_y,               \
                  #T, (int)CK_FLOATING_DIG, _ck_t);                      \
    }                                                                    \
  } while (0)


/* Floating point relative tolerance comparison for arrays. Based on
 * ck_assert_floating_absdiff_tol
 * OP can have values: <; >=.

 * Should print an failure message like
 * Assertion '|u[0] - u_exp[0]|/min(|u[0]| - |u_exp[0]|) < 1e-14' failed:
 * u_exp[0]  == 1.3500, u[0] == 10000,      min(|u_exp[0]|, |u[0]|)) == 1.3500,  TOL == 1e-14
 * u[0]/min(|u|, |u_exp|) = 7407.4074074074069,u_exp[0]/min(|u|, |u_exp|) = 1,
 */
#define _ck_assert_floating_array_absdiff_op_reltol(N, X, Y, OP, T, TP)           \
  do {                                                                            \
    TP _ck_rel_num=1;                                                             \
                                                                                  \
    for (int _ck_i_=0; _ck_i_<N;_ck_i_++){                                        \
      TP _ck_x = (X[_ck_i_]);                                                     \
      TP _ck_y = (Y[_ck_i_]);                                                     \
      _ck_rel_num = ck_min(fabsl(_ck_x), fabsl(_ck_y));                           \
      TP _ck_t = (T);                                                             \
      ck_assert_msg(fabsl(_ck_y - _ck_x)/_ck_rel_num OP _ck_t,                    \
                    "\n     Assertion '|%s[%d] - %s[%d]|"                         \
                    "/min(|%s[%d]| - |%s[%d]|) %s %s' failed: \n"                 \
                    "     %s[%d]  == %.*g,"                                       \
                    " %s[%d] == %.*g,\n "                                         \
                    "    min(|%s[%d]|, |%s[%d]|)) == %.*g,  "                     \
                    "TOL == %.*g\n"                                               \
                    "     %s[%d]/min(|%s|, |%s|) = %.*g, "                        \
                    "%s[%d]/min(|%s|, |%s|) = %.*g\n",                            \
                    #Y,_ck_i_, #X,_ck_i_, #Y,_ck_i_, #X,_ck_i_, #OP, #T,          \
                    #X, _ck_i_, (int)CK_FLOATING_DIG, _ck_x,                      \
                    #Y, _ck_i_, (int)CK_FLOATING_DIG, _ck_y,                      \
                    #X, _ck_i_, #Y, _ck_i_, (int)CK_FLOATING_DIG, _ck_rel_num,    \
                    (int)CK_FLOATING_DIG, _ck_t,                                  \
                    #Y, _ck_i_, #Y, #X, (int)CK_FLOATING_DIG, _ck_y/_ck_rel_num,  \
                    #X, _ck_i_, #Y, #X, (int)CK_FLOATING_DIG, _ck_x/_ck_rel_num); \
    }                                                                             \
  } while (0)



/*
 * Check two double precision floating point arrays to determine
 * if all elements satisfy the relative tolerance condition
 * if |X[i] - Y[i]|/min(|X[i}, |Y[i]|) < TOL
 * for i=0,...N-1
 *
 * If |X[i] - Y[i]|/min(|X[i}, |Y[i]|) >= T for any i, the test fails.
 *
 * @param N Number of elements in array.
 * @param X floating point number (double)
 * @param Y floating point number (double) to compare against X
 * @param T tolerance (double)
 *
 * @note If the check fails, the remaining of the test is aborted
 *
 * @since Not part of standard libcheck
 */

#define ck_assert_double_array_eq_reltol(N, X, Y, T)  \
  _ck_assert_floating_array_absdiff_op_reltol(N, X, Y, <, T, double)
/*
 * Check two double precision floating point array to determine if X[i]/nrm≈Y[i]/nrm
 * with specified tolerance, for i=0,...N-1
 *
 */



/*
 * Check two double precision floating point arrays to determine if |X[i] - Y[i]| < TOL
 * for i=0,...N-1
 *
 * If |X[i] - Y[i]| >= T for any i, the test fails.
 *
 * @param N Number of elements in array.
 * @param X floating point number (double)
 * @param Y floating point number (double) to compare against X
 * @param T tolerance (double)
 *
 * @note If the check fails, the remaining of the test is aborted
 *
 * @since Not part of standard libcheck
 */
#define ck_assert_double_array_eq_tol(N, X, Y, T)  \
  _ck_assert_floating_array_absdiff_op_tol(N, X, Y, <, T, double, "")

/*
 * Check two double precision floating point arrays to determine if X[i] == Y[i]
 * for i=0,...N-1
 *
 * If X[i] != Y[i] for any i, the test fails.
 *
 * @param N Number of elements in array.
 * @param X floating point number (double)
 * @param Y floating point number (double) to compare against X
 * @param T tolerance (double)
 *
 * @note If the check fails, the remaining of the test is aborted
 *
 * @since Not part of standard libcheck
 */
#define ck_assert_double_array_eq(N, X, Y)                       \
  _ck_assert_floating_array_absdiff_op_tol(N, X, Y, <=, 0, double, "")
/* using <= with T=0 prevents warnings about equality comparison of double.*/


#endif

