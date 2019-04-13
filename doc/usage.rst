=====
Usage
=====

This page describes writing `c` programs which call `L1c`. For usage of the Python or Matlab bindings, see :ref:`bindings`.

The following is incomplete and only describes the main library interface to the `l1qc` solver.

As a user, the primary function you need to worry about is

.. code-block:: c

   /*l1qc_newton.h */
   LBResult l1qc_newton(l1c_int N, double *x, l1c_int M, double *b,
                NewtParams params, L1cAxFuns Ax_funs);


* `int N`. The length of `x` and `u`.
* `double *x`. On entry, this should be an array of doubles length N, allocated on a 64-byte boundary (see below). On exit, x contains the result of the optimization.
* `double *u` On entry, this should contain an array with length N. On exit, it will contain the auxilary u (See above about the conversion from an l1 optimization to a linear program).
* `int M`. The length b.
* `double *b`. On entry, contains the 'measured data' (see above). In general, we expect M <N.
* `NewtParams params` is a struct containing parameters (e.g., tolerances and iteration number bounds). Will be described fully below.
* `L1cAxFuns ax_funs` is a struct containing pointers to the functions which perform the transformations.


*Important*: The array inputs of doubles (`*x, *u, *b`) to `l1qc_newton` must be aligned on a 64-byte boundary, otherwise segfaults may occur. To faciliate this, you may use the functions 

.. code-block:: c

   /*l1c_common.h */
   void* malloc_double(N);
   void* free_double(N);

The function :c:func:`malloc_double` will allocate memory for `N` doubles, aligned on a 64-byte boundary and `free_double` will free it.


The data structures are defined as

.. code-block:: c

   //l1qc_newton.h
   typedef struct LBResult{
                double l1;                // Final value of functional, ||x||_1
                int    total_newton_iter; // Total number of newton iterations.
                int    status;            // 0 if completed with no errors, 1 otherwise

                }LBResult;

   typedef struct NewtParams{
                double epsilon;
                double tau;
                double mu;
                double newton_tol;
                int newton_max_iter;
                int lbiter;
                double lbtol;
                int verbose;
                CgParams cg_params;

   }NewtParams;

   typedef struct L1cAxFuns {
                void(*Ax)(double *x, double *y);
                void(*Aty)(double *y, double *x);
                void(*AtAx)(double *x, double *z);
                }L1cAxFuns;


The struct `L1cAxFuns` contains pointers to your user-defined functions which compute latexmath:[Ax] and latexmath:[A^{T}y] For an example, see the mex-interface file `l1qc_mex.c` (in `interfaces/`) and `dct1.c` or `dct2.c`. Note that although the mex interface looks long and complicated, almost all of this is boiler-plate parsing of Matlab's input to the function. The amount of code to modify for a different set of transform functions is only a few lines.
