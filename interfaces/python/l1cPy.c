#include "config.h"
#include <stdlib.h>

#include "l1c.h"
#include "nesta.h"
#include <stdio.h>
#include <Python.h>
// Without the following, the numpy header generates #warnings
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/ndarraytypes.h>



/* Function declartion */
static PyObject *_l1qc_dct(PyObject *self, PyObject *args, PyObject *kw);
static PyObject *_breg_anistropic_TV(PyObject *self, PyObject *args, PyObject *kw);
static PyObject *_nesta_dctTV(PyObject *self, PyObject *args, PyObject *kw);

PyMODINIT_FUNC PyInit__l1cPy_module(void);

/* Doc-strings */
static char module_docstring[] =
  "This module provides an interface for solving a l1 optimization problems in c.";
static char l1qc_dct_docstring[] =  "Minimize ||x||_1 s.t. ||Ax-b|| < epsilon";

static char nesta_dctTV_docstring[] =
  "Minimize ||U^TF||_1 + alp_h||\\nabla_xF|| + alp_v||\\nabla_v F||s.t. ||Ax-b|| < sigma";

static char breg_anistropic_TV_docstring[] =
  "Given an n by m image f, solves the anistropic TV denoising problem \n"
  "f_denoised = breg_anistropic_TV(f, mu=10, tol=0.001, max_iter=1000, max_jac_iter=1)\n"
  " \n"
  "min ||\\nabla_x u||_1 + ||\\nabla_y||_1 + 0.5\\mu ||u - f||_2 \n"
   "u \n"
  " \n"
  "using Bregman Splitting. This algorithm is developed in the paper \n"
  "\"The Split Bregman Method for L1-Regularized Problems,\" by Tom Goldstein and Stanley Osher. ";


/* Specify members of the module */
static PyMethodDef _l1cPy_module_methods[] = {
  {"l1qc_dct",
   (PyCFunction)_l1qc_dct,
   METH_KEYWORDS|METH_VARARGS,
   l1qc_dct_docstring},

  {"breg_anistropic_TV",
   (PyCFunction)_breg_anistropic_TV,
   METH_KEYWORDS|METH_VARARGS,
   breg_anistropic_TV_docstring},
   {"nesta_dctTV",
   (PyCFunction)_nesta_dctTV,
   METH_KEYWORDS|METH_VARARGS,
    nesta_dctTV_docstring},

  {NULL, NULL, 0, NULL}
};


static struct PyModuleDef mod_l1cPy = {
    PyModuleDef_HEAD_INIT,
    "_l1cPy_module",          /* name of module */
    module_docstring,         /* module documentation */
    -1,                       /* size of per-interpreter state of the module,
                               or -1 if the module keeps state in global variables. */
    _l1cPy_module_methods,    /* m_methods */
    /* Zero the rest to avoid warnings. */
    NULL,                     /* m_slots*/
    0,                        /* traverseproc, int in object.h*/
    0,                        /* m_clear*/
    NULL,                     /* m_free*/

  };

PyMODINIT_FUNC PyInit__l1cPy_module(void)
{
  import_array();
  return PyModule_Create(&mod_l1cPy);
}


static PyObject *
_l1qc_dct(PyObject *self, PyObject *args, PyObject *kw){
  (void)self;
  double *x_ours=NULL, *b=NULL;

  PyObject *b_obj=NULL, *pix_idx_obj=NULL, *x_out_npA=NULL;
  PyObject *lb_res_obj=NULL, *ret_val=NULL;
  PyArrayObject *b_npA=NULL;
  PyArrayObject *pix_idx_npA=NULL;

  int i=0, status=0, mrow=0, mcol=0, mtot=0, n=0, n_=0;
  int *pix_idx=NULL;

  l1c_LBResult lb_res = {.status = 0, .total_newton_iter = 0, .l1=INFINITY};

  l1c_L1qcOpts opts = {.epsilon=.01, .mu = 10,
                       .lbtol = 1e-3, .newton_tol = 1e-3,
                       .newton_max_iter = 50, .verbose = 0,
                       .l1_tol = 1e-5, .lbiter = 0,
                       .cg_tol = 1e-8, .cg_maxiter = 200,
                       .cg_verbose=0, .warm_start_cg=0, .dct_mode=dct1};

  char *keywords[] = {"", "", "", "",
                      "epsilon", "mu",
                      "lbtol", "newton_tol",
                      "newton_max_iter", "verbose",
                      "l1_tol", "cgtol",
                      "cgmaxiter", "dct_mode",
                      NULL};

  /* Parse the input tuple O=object, d=double, i=int*/
  if (!PyArg_ParseTupleAndKeywords(args, kw, "iiOO|ddddiiddii", keywords,
                                   &mrow, &mcol, &b_obj, &pix_idx_obj,
                                   &opts.epsilon, &opts.mu,
                                   &opts.lbtol, &opts.newton_tol,
                                   &opts.newton_max_iter, &opts.verbose,
                                   &opts.l1_tol, &opts.cg_tol,
                                   &opts.cg_maxiter, &opts.dct_mode)){
    PyErr_SetString(PyExc_TypeError, "Parsing input arguments failed.");
    return NULL;
  }

  if (opts.dct_mode != dct1 && opts.dct_mode != dct2) {
    PyErr_SetString(PyExc_ValueError, "Must have dct_mode = 1 or 2");
    return NULL;
  }

  /* Interpret the input objects as numpy arrays.
     N.B: NPY_ARRAY_IN_ARRAY = PY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED
   */
  b_npA =(PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  pix_idx_npA =(PyArrayObject*)PyArray_FROM_OTF(pix_idx_obj, NPY_INT, NPY_ARRAY_IN_ARRAY);


  /* If that didn't work, throw an exception. */
  if ( b_npA == NULL || pix_idx_npA == NULL ){
    PyErr_SetString(PyExc_RuntimeError, "Failed to initialize arrays");
    goto fail;
  }

  /*--------------- Check number of Dimensions ----------------- */
  if ( (PyArray_NDIM(b_npA) != 1) || (PyArray_NDIM(pix_idx_npA) != 1)){
    PyErr_SetString(PyExc_IndexError, "b and pix_idx must have exactly 1 dimensions.");
    goto fail;
  }

  /*--------------- Check Sizes ----------------- */
  mtot = mrow*mcol;
  n = PyArray_DIM(b_npA, 0);
  n_ = PyArray_DIM(pix_idx_npA, 0);

  if ( mtot < n || n != n_ ){
    PyErr_SetString(PyExc_ValueError, "Must have len(x) > len(b) and len(b) == len(pix_idx).");
    goto fail;
  }

  /* Get pointers to the data as C-types. */
  b  = (double*)PyArray_DATA(b_npA);
  pix_idx = (int*)PyArray_DATA(pix_idx_npA);

  /* Allocate memory for M*xk=f */
  x_ours = l1c_malloc_double(mtot);
  if (!x_ours){
    PyErr_SetString(PyExc_MemoryError, "Failed to allocation memory");
    goto fail;
  }

  for(i=0; i<mtot; i++){
    x_ours[i] = 0;
  }

  status = l1qc_dct(mrow, mcol, x_ours, n, b, pix_idx, opts, &lb_res);
  if (status != L1C_SUCCESS){
    PyErr_SetString(PyExc_RuntimeError, "l1qc_dct failed.");
    l1c_free_double(x_ours);
    goto fail;
  }
  /* Build the output tuple */
  npy_intp dims[] = {mrow, mcol};
  x_out_npA = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, x_ours);
  lb_res_obj = Py_BuildValue("{s:d,s:i,s:i,s:i}",
                                       "l1", lb_res.l1,
                                       "total_newton_iter", lb_res.total_newton_iter,
                                       "total_cg_iter", lb_res.total_cg_iter,
                                       "status", lb_res.status|status);

  /* Clean up. */
  Py_DECREF(b_npA);
  Py_DECREF(pix_idx_npA);

  ret_val = Py_BuildValue("OO", x_out_npA, lb_res_obj);
  return ret_val;


  /* If we failed, clean up more things */
 fail:
  Py_XDECREF(b_npA);
  Py_XDECREF(pix_idx_npA);

  return NULL;

}



static PyObject *
_breg_anistropic_TV(PyObject *self, PyObject *args, PyObject *kw){
  (void)self;

  PyObject *f_obj=NULL, *uk_out_npA=NULL;
  PyObject *ret_val=NULL;
  PyArrayObject *f_npA=NULL;
  double *uk=NULL, *f=NULL, *f_ours=NULL;
  int i=0, status=0, n=0, m=0;

  double tol=0.001, mu=10;
  l1c_int max_iter = 1000;
  l1c_int max_jac_iter = 1;

  char *keywords[] = {"",
                      "mu", "tol",
                      "max_iter", "max_jac_iter",
                      NULL};

  /* Parse the input tuple O=object, d=double, i=int*/
  if (!PyArg_ParseTupleAndKeywords(args, kw, "O|ddii", keywords,
                                   &f_obj,
                                   &mu, &tol, &max_iter,&max_jac_iter)){
    PyErr_SetString(PyExc_TypeError, "Parsing input arguments failed.");
    return NULL;
  }
  /* Interpret the input objects as numpy arrays.
     N.B: NPY_ARRAY_IN_ARRAY = PY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED
   */
  f_npA =(PyArrayObject*)PyArray_FROM_OTF(f_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

  /* If that didn't work, throw an exception. */
  if ( f_npA == NULL){
    PyErr_SetString(PyExc_RuntimeError, "Failed to initialize arrays");
    goto fail1;
  }

  /*--------------- Check number of Dimensions ----------------- */
  if (PyArray_NDIM(f_npA) != 2){
    PyErr_SetString(PyExc_IndexError, "f must have exactly 2 dimensions.");
    goto fail1;
  }

  /*--------------- Check Sizes ----------------- */
  n = PyArray_DIM(f_npA, 0);
  m = PyArray_DIM(f_npA, 1);
  if ( n <3 || m <3 ){
    PyErr_SetString(PyExc_ValueError, "Must have both dimensions of f greater than 2.");
    goto fail1;
  }

  /* Get pointers to the data as C-types. */
  f  = (double*)PyArray_DATA(f_npA);

  /* Allocate memory for M*xk=f */
  f_ours = l1c_malloc_double(n*m);
  uk = l1c_malloc_double(n*m);
  if (!f_ours || !uk){
    PyErr_SetString(PyExc_MemoryError, "Failed to allocation memory");
    goto fail1;
  }

  for(i=0; i<n*m; i++){
    f_ours[i] = f[i];
  }

  status = l1c_breg_anistropic_TV(n, m, uk, f_ours, mu, tol, max_iter, max_jac_iter);
  if (status){
    PyErr_SetString(PyExc_MemoryError, "Failed to allocation memory");
    goto fail2;
  }
  /* Build the output tuple */
  npy_intp dims[] = {n, m};
  uk_out_npA = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, uk);

  /* Clean up. */
  Py_DECREF(f_npA);
  l1c_free_double(f_ours);

  ret_val = Py_BuildValue("O", uk_out_npA);

  return ret_val;

  /* If we failed, clean up more things. Note the fall through. */
 fail2:
  l1c_free_double(f_ours);
  l1c_free_double(uk);
 fail1:
  Py_XDECREF(f_npA);

  return NULL;

}



static PyObject *
_nesta_dctTV(PyObject *self, PyObject *args, PyObject *kw){
  (void)self;
  double *x_ours=NULL, *b=NULL;

  PyObject *b_obj=NULL, *pix_idx_obj=NULL, *x_out_npA=NULL;

  PyArrayObject *b_npA=NULL;
  PyArrayObject *pix_idx_npA=NULL;

  int status=0, mrow=0, mcol=0, mtot=0, n=0, n_=0;
  int *pix_idx=NULL;
  double alpha_v=0, alpha_h=0;
  l1c_NestaOpts opts = {.sigma = 0.01, .mu = 10,
                        .verbose = 5, .tol = 1e-5,
                        .n_continue=5};

  DctMode dct_mode = dct1;
  BpMode bp_mode = analysis;
  // int _dct_mode, _bp_mode;
  char *keywords[] = {"", "", "", "",
                      "alpha_v", "alpha_h",
                      "sigma", "mu", "tol", "n_continue", "verbose",
                      "dct_mode", "bp_mode",
                      NULL};

  /* Parse the input tuple O=object, d=double, i=int*/
  if (!PyArg_ParseTupleAndKeywords(args, kw, "iiOO|dddddiiii", keywords,
                                   &mrow, &mcol, &b_obj, &pix_idx_obj,
                                   &alpha_v, &alpha_h,
                                   &opts.sigma, &opts.mu,
                                   &opts.tol, &opts.n_continue, &opts.verbose,
                                   &dct_mode, &bp_mode)){

    PyErr_SetString(PyExc_TypeError, "Parsing input arguments failed.");
    return NULL;
  }
  if (bp_mode != analysis && bp_mode != synthesis) {
    PyErr_SetString(PyExc_ValueError, "Must have bp_mode = 1 or 2");
    return NULL;
  }

  if (dct_mode != dct1 && dct_mode != dct2) {
    PyErr_SetString(PyExc_ValueError, "Must have dct_mode = 1 or 2");
    return NULL;
  }

  /* Interpret the input objects as numpy arrays.
     N.B: NPY_ARRAY_IN_ARRAY = PY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED */

  b_npA =(PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  pix_idx_npA =(PyArrayObject*)PyArray_FROM_OTF(pix_idx_obj, NPY_INT, NPY_ARRAY_IN_ARRAY);

  /* If that didn't work, throw an exception. */
  if ( b_npA == NULL || pix_idx_npA == NULL ){
    PyErr_SetString(PyExc_RuntimeError,
        "Failed to initialize arrays. Is pix_idx an integer array?");
    goto fail;
  }


  /*--------------- Check number of Dimensions ----------------- */
  if ( (PyArray_NDIM(b_npA) != 1) || (PyArray_NDIM(pix_idx_npA) != 1)){
    PyErr_SetString(PyExc_IndexError,
                    "b and pix_idx must have exactly 1 dimensions.");
    goto fail;
  }


  /*--------------- Check Sizes ----------------- */
  mtot = mrow*mcol;
  n = PyArray_DIM(b_npA, 0);
  n_ = PyArray_DIM(pix_idx_npA, 0);

  if ( mtot < n || n != n_ ){
    PyErr_SetString(PyExc_ValueError, "Must have len(x) > len(b) and len(b) == len(pix_idx).");
    goto fail;
  }

  /* Get pointers to the data as C-types. */
  b  = (double*)PyArray_DATA(b_npA);
  pix_idx = (int*)PyArray_DATA(pix_idx_npA);

  /* Allocate memory for M*xk=f */
  x_ours = l1c_calloc_double(mtot);
  if (!x_ours){
    PyErr_SetString(PyExc_MemoryError, "Failed to allocation memory");
    goto fail;
  }


  l1c_AxFuns ax_funs;
  if (l1c_setup_dctTV_transforms(n, mrow, mcol, alpha_v, alpha_h, dct_mode,
                                 bp_mode, pix_idx, &ax_funs)) {
    PyErr_SetString(PyExc_RuntimeError, "Failed to initialize DCT.");
    goto fail;
  }

  status = l1c_nesta(mtot, x_ours, n, b, ax_funs, opts);

  /* Build the output tuple */
  npy_intp dims[] = {mrow, mcol};
  x_out_npA = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, x_ours);

  /* Clean up. */
  Py_DECREF(b_npA);
  Py_DECREF(pix_idx_npA);

  return Py_BuildValue("Oi", x_out_npA, status);

  /* If we failed, clean up more things */
 fail:
  Py_XDECREF(b_npA);
  Py_XDECREF(pix_idx_npA);

  return NULL;
}

