#include "config.h"
#include <stdlib.h>
#include <cblas.h>

#include "l1c.h"
#include "nesta.h"
#include <stdio.h>
#include <Python.h>
// Without the following, the numpy header generates #warnings
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/ndarraytypes.h>

/* Helper functions*/
static PyObject *PyArray_SimpleDoubleCopy(int nd, npy_intp *dims, double *data);
static int is_bpmode_valid(BpMode bpmode);
static int is_dctmode_valid(DctMode dctmode);

/* Exported Function declartion */
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
  PyObject *lb_res_obj=NULL;
  PyArrayObject *b_npA=NULL;
  PyArrayObject *pix_idx_npA=NULL;

  int status=0, mrow=0, mcol=0, mtot=0, n=0, n_=0;
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

  if (!is_dctmode_valid(opts.dct_mode))
    return NULL;

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
  x_ours = l1c_calloc_double(mtot);
  if (!x_ours){
    PyErr_SetString(PyExc_MemoryError, "Failed to allocation memory");
    goto fail;
  }


  status = l1qc_dct(mrow, mcol, x_ours, n, b, pix_idx, opts, &lb_res);
  if (status != L1C_SUCCESS){
    PyErr_SetString(PyExc_RuntimeError, "l1qc_dct failed.");
    goto fail;
  }
  /* Build the output tuple */
  npy_intp dims[] = {mrow, mcol};
  // x_out_npA = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, x_ours);
  x_out_npA = PyArray_SimpleDoubleCopy(2, dims, x_ours);
  if (!x_out_npA)
    goto fail;

  lb_res_obj = Py_BuildValue("{s:d,s:i,s:i,s:i}",
                                       "l1", lb_res.l1,
                                       "total_newton_iter", lb_res.total_newton_iter,
                                       "total_cg_iter", lb_res.total_cg_iter,
                                       "status", lb_res.status|status);
  if(!lb_res_obj)
    goto fail2;

  /* Clean up. */
  Py_DECREF(b_npA);
  Py_DECREF(pix_idx_npA);
  l1c_free_double(x_ours);

  return Py_BuildValue("OO", x_out_npA, lb_res_obj);

  /* If we failed, clean up more things */
 fail2:
   Py_DECREF(x_out_npA);
 fail:
  Py_XDECREF(b_npA);
  Py_XDECREF(pix_idx_npA);
  l1c_free_double(x_ours);

  return NULL;

}



static PyObject *
_breg_anistropic_TV(PyObject *self, PyObject *args, PyObject *kw){
  (void)self;

  PyObject *f_obj=NULL, *uk_out_npA=NULL;

  PyArrayObject *f_npA=NULL;
  double *uk=NULL, *f=NULL, *f_ours=NULL;
  int status=0, n=0, m=0;

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
  if (!f_ours || !uk || !f){
    PyErr_SetString(PyExc_MemoryError, "Failed to allocation memory");
    goto fail1;
  }

  cblas_dcopy(n*m, f, 1, f_ours, 1);

  status = l1c_breg_anistropic_TV(n, m, uk, f_ours, mu, tol, max_iter, max_jac_iter);
  if (status){
    PyErr_SetString(PyExc_MemoryError, "Failed to allocation memory");
    goto fail2;
  }
  /* Build the output tuple */
  npy_intp dims[] = {n, m};
  uk_out_npA = PyArray_SimpleDoubleCopy(2, dims, uk);
  if(!uk_out_npA)
    goto fail2;

  /* Clean up. */
  Py_DECREF(f_npA);
  l1c_free_double(f_ours);
  l1c_free_double(uk);

  return Py_BuildValue("O", uk_out_npA);

  /* If we failed, clean up more things. Note the fall through. */
 fail2:
  l1c_free_double(f_ours);
  l1c_free_double(uk);
 fail1:
  Py_XDECREF(f_npA);

  return NULL;

}


static PyObject *_nesta_dctTV(PyObject *self, PyObject *args,
                                                PyObject *kw) {
  (void)self;
  double *x_ours=NULL, *b=NULL, *b_ours=NULL;

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
  l1c_AxFuns ax_funs;

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

  if (!is_bpmode_valid(bp_mode))
    return NULL;

  if (!is_dctmode_valid(dct_mode))
    return NULL;


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
  b_ours = l1c_malloc_double(n);
  if (!x_ours || !b_ours){
    PyErr_SetString(PyExc_MemoryError, "Failed to allocation memory");
    goto fail;
  }

  cblas_dcopy(n, b, 1, b_ours, 1);


  if (l1c_setup_dctTV_transforms(n, mrow, mcol, alpha_v, alpha_h, dct_mode,
                                 bp_mode, pix_idx, &ax_funs)) {
    PyErr_SetString(PyExc_RuntimeError, "Failed to initialize DCT.");
    goto fail;
  }

  status = l1c_nesta(mtot, x_ours, n, b_ours, ax_funs, opts);

  /* Build the output tuple */
  npy_intp dims[] = {mrow, mcol};

  x_out_npA = PyArray_SimpleDoubleCopy(2, dims, x_ours);
  if (!x_out_npA) goto fail;

  /* Clean up. */
  Py_DECREF(b_npA);
  Py_DECREF(pix_idx_npA);
  l1c_free_double(b_ours);
  l1c_free_double(x_ours);

  return Py_BuildValue("Oi", x_out_npA, status);

  /* If we failed, clean up more things */
 fail:
  l1c_free_double(x_ours);
  l1c_free_double(b_ours);
  Py_XDECREF(b_npA);
  Py_XDECREF(pix_idx_npA);

  return NULL;
}

static int is_bpmode_valid(BpMode bp_mode){
if (bp_mode != analysis && bp_mode != synthesis) {
  PyErr_SetString(PyExc_ValueError, "Must have bp_mode = 1 or 2");
  return 0;
 }

 return 1;
}

static int is_dctmode_valid(DctMode dct_mode){
  if (dct_mode != dct1 && dct_mode != dct2) {
    PyErr_SetString(PyExc_ValueError, "Must have dct_mode = 1 or 2");
    return 0;
  }

  return 1;
}

/*
  When building an output array from data we allocated, we cant free our own
  memory, because numpy just uses a reference to it, but cant decref properly.
  This creates a memory leak. So the point of this function is to create an
  output array, and copy our data into it, so we can free our own memory.


  A better idea is:
  http://blog.enthought.com/python/numpy-arrays-with-pre-allocated-memory/

  which shows how to register a deallocator with numpy.

  See also,
  https://docs.scipy.org/doc/numpy/reference/c-api.array.html#c.PyArray_New
 */
static PyObject *PyArray_SimpleDoubleCopy(int nd, npy_intp *dims, double *data) {
  int N = 1;
  PyObject *py_arr_obj = PyArray_SimpleNew(nd, dims, NPY_DOUBLE);
  if (!py_arr_obj) {
    return NULL;
  }

  PyArrayObject *arr_npA = (PyArrayObject *)PyArray_FROM_OTF(
      py_arr_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  if (!arr_npA) {
    goto fail1;
  }

  double *their_data = (double *)PyArray_DATA(arr_npA);
  if (!their_data) {
    goto fail2;
  }
  for (int i = 0; i < nd; i++) {
    N *= dims[i];
  }

  cblas_dcopy(N, data, 1, their_data, 1);

  Py_XDECREF(arr_npA);

  return py_arr_obj;

fail2:
  Py_XDECREF(arr_npA);
fail1:
  Py_XDECREF(py_arr_obj);
  return NULL;
}
