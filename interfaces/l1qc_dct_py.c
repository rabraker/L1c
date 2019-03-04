#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/ndarraytypes.h>

#include "fgm.h"
#include "utils.h"
#include "cblas.h"
#include <stdlib.h>
#include <stdio.h>


/* Function declartion */
static PyObject *l1qc(PyObject *self, PyObject *args);
PyMODINIT_FUNC PyInit_fgmPy(void);


/* Doc-strings */
static char module_docstring[] =
  "This module provides an interface for solving a l1 optimization problems in c.";
static char l1qc_docstring[] =
  "Minimize ||x||_1 s.t. ||Ax-b|| < epsilon";



/* Specify members of the module */
static PyMethodDef l1cPy_methods[] = {
  {"l1qc_dct", l1qc_dct, METH_VARARGS, l1qc_docstring},
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


static PyObject *l1qc_dct(PyObject *self, PyObject *args)
{

  PyObject *x_obj, *b_obj, *pix_idx_obj;
  double umax, beta;
  int maxIter;

  /* Parse the input tuple O=object, d=double, i=int*/
    if (!PyArg_ParseTuple(args, "OOOOOddi", &IH_obj, &M_obj, &xk_obj, &zi_obj, &yi_obj,
                          &beta, &umax, &maxIter)){
    printf("Fail 000\n");
    return NULL;
  }
  /* Interpret the input objects as numpy arrays. */
  PyArrayObject *x_array =(PyArrayObject*)PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  PyArrayObject *b_array =(PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  PyArrayObject *pix_idx_array =(PyArrayObject*)PyArray_FROM_OTF(pix_idx_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);


  /* If that didn't work, throw an exception. */
  if (x_array == NULL  || b_array == NULL || pix_idx_array == NULL ||){
    PyErr_SetString(PyExc_RuntimeError, "Failed to initialize arrays");
    goto fail;
  }

  /*--------------- Check number of Dimensions ----------------- */
  if ( (PyArray_NDIM(x_array) != 1) || (PyArray_NDIM(b_array) != 1)
       ||(PyArray_NDIM(pix_idx_array) != 1)){
    PyErr_SetString(PyExc_IndexError, "IH and M must have exactly 1 dimensions.");
    goto fail;
  }

  /*--------------- Check Sizes ----------------- */
  int N = PyArray_DIM(x_array, 0);
  int M = PyArray_DIM(b_array, 0);
  int M_ = PyArray_DIM(pix_idx_array, 1);

  if ( N < M || M != M_ ){
    PyErr_SetString(PyExc_IndexError, "Must have len(x) > len(b) \
and len(b) == len(pix_idx).");
    goto fail;
  }



  /* Get pointers to the data as C-types. */
  double *x_theirs = (double*)PyArray_DATA(x_array);
  double *b  = (double*)PyArray_DATA(b_array);
  int *pix_idx = (double*)PyArray_DATA(pix_idx_array);

  /* Allocate memory for M*xk=f */
  double *x_ours;
  x_ours = malloc_double(N)
  if (!x_ours){
    perror("Memory Allocation failure\n");
    free(x_ours);
    PyErr_SetString(PyExc_MemoryError, "Failed to allocation memory");
  }
  /* Pointer to a function that computes Ax=f. Python passes us the matrix in
     RowMajor order. Signature: dgemv_RowOrder(M,m_M, n_M, xk, f)
  */
  void (*AX_func)(double *, int, int, double *, double *) = dgemv_RowOrder;
  AX_func(M, m_M, n_M, xk, f);

  if (solve_fgm(IH, m_IH, f, beta, maxIter, umax, zi, yi, AX_func) ){goto fail;}


  /* Clean up. */
  Py_DECREF(IH_array);
  Py_DECREF(M_array);
  Py_DECREF(xk_array);
  free(f);

  /* Build the output tuple */
  PyObject *ret = Py_BuildValue("OO", zi_array, yi_array);
  return ret;
  /* If we failed, clean up more things */
 fail:
  Py_XDECREF(IH_array);
  Py_XDECREF(M_array);
  Py_XDECREF(xk_array);
  Py_XDECREF(zi_array);
  Py_XDECREF(yi_array);
  printf("INSIDE\n");
  return NULL;

}

