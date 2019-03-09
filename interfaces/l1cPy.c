#include "config.h"
#include "l1c_common.h"
#include <stdlib.h>
#include <stdio.h>

#include "Python.h"
// Without the following, the numpy header generates #warnings
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/ndarraytypes.h>

#include "l1qc_newton.h"
#include "libl1qc_dct.h"


/* Function declartion */
static PyObject *l1qc_dct_py(PyObject *self, PyObject *args, PyObject *kw);
PyMODINIT_FUNC PyInit_l1cPy(void);


/* Doc-strings */
static char module_docstring[] =
  "This module provides an interface for solving a l1 optimization problems in c.";
static char l1qc_docstring[] =  "Minimize ||x||_1 s.t. ||Ax-b|| < epsilon";



/* Specify members of the module */
static PyMethodDef l1cPy_methods[] = {
  {"l1qc_dct_py", (PyCFunction)l1qc_dct_py, METH_KEYWORDS|METH_VARARGS, l1qc_docstring},
  {NULL, NULL, 0, NULL}
};


static struct PyModuleDef mod_l1cPy = {
    PyModuleDef_HEAD_INIT,
    "l1cPy",          /* name of module */
    module_docstring, /* module documentation */
    -1,               /* size of per-interpreter state of the module,
                         or -1 if the module keeps state in global variables. */
    l1cPy_methods
  };

PyMODINIT_FUNC PyInit_l1cPy(void)
{
  import_array();
  return PyModule_Create(&mod_l1cPy);
}


static PyObject *
l1qc_dct_py(PyObject *self, PyObject *args, PyObject *kw){

  double *x_ours=NULL, *b=NULL;

  PyObject *x_obj=NULL, *b_obj=NULL, *pix_idx_obj=NULL, *x_out_npA=NULL;
  PyObject *lb_res_obj=NULL, *ret_val=NULL;
  PyArrayObject *b_npA=NULL;
  PyArrayObject *pix_idx_npA=NULL;

  int i=0, status=0, n=0, m=0, N=0, M=0, M_=0;
  int *pix_idx=NULL;

  LBResult lb_res = {.status = 0, .total_newton_iter = 0, .l1=INFINITY};

  L1qcDctOpts opts = {.epsilon=.01, .mu = 10,
                      .lbtol = 1e-3, .newton_tol = 1e-3,
                      .newton_max_iter = 50, .verbose = 0,
                      .l1_tol = 1e-5, .lbiter = 0,
                      .cgtol = 1e-8, .cgmaxiter = 200,
                      .warm_start_cg=0};

  char *keywords[] = {"", "", "", "",
                      "epsilon", "mu",
                      "lbtol", "newton_tol",
                      "newton_max_iter", "verbose",
                      "l1_tol", "cgtol",
                      "cgmaxiter",
                      NULL};

  /* Parse the input tuple O=object, d=double, i=int*/
  if (!PyArg_ParseTupleAndKeywords(args, kw, "iiOO|ddddiiddi", keywords,
                                   &n, &m, &b_obj, &pix_idx_obj,
                                   &opts.epsilon, &opts.mu,
                                   &opts.lbtol, &opts.newton_tol,
                                   &opts.newton_max_iter, &opts.verbose,
                                   &opts.l1_tol, &opts.cgtol,
                                   &opts.cgmaxiter)){
    fprintf(stderr, "Parsing input arguments failed.\n");
    return NULL;
  }

  /* Interpret the input objects as numpy arrays. */
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
  N = n*m;
  M = PyArray_DIM(b_npA, 0);
  M_ = PyArray_DIM(pix_idx_npA, 0);

  if ( N < M || M != M_ ){
    PyErr_SetString(PyExc_ValueError, "Must have len(x) > len(b) and len(b) == len(pix_idx).");
    goto fail;
  }

  /* Get pointers to the data as C-types. */
  b  = (double*)PyArray_DATA(b_npA);
  pix_idx = (int*)PyArray_DATA(pix_idx_npA);

  /* Allocate memory for M*xk=f */
  x_ours = malloc_double(N);
  if (!x_ours){
    fprintf(stderr, "Memory Allocation failure\n");
    PyErr_SetString(PyExc_MemoryError, "Failed to allocation memory");
  }

  for(i=0; i<N; i++){
    x_ours[i] = 0;
  }


  status = l1qc_dct(N, 1, x_ours, M, b, pix_idx, opts, &lb_res);


  /* Clean up. */
  Py_DECREF(b_npA);
  Py_DECREF(pix_idx_npA);

  /* Build the output tuple */
  npy_intp dims[] = {N};
  x_out_npA = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, x_ours);
  lb_res_obj = Py_BuildValue("{s:d,s:i,s:i,s:i}",
                                       "l1", lb_res.l1,
                                       "total_newton_iter", lb_res.total_newton_iter,
                                       "total_cg_iter", lb_res.total_cg_iter,
                                       "status", lb_res.status);

  ret_val = Py_BuildValue("OO", x_out_npA, lb_res_obj);
  return ret_val;


  /* If we failed, clean up more things */
 fail:
  Py_DECREF(b_npA);
  Py_DECREF(pix_idx_npA);

  // Py_XDECREF(x_ours);
  printf("INSIDE\n");
  return NULL;

}

