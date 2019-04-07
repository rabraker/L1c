import numpy as np
from ctypes import Structure
from ctypes import c_int
from ctypes import c_double
from ctypes import POINTER, byref
import time


def BRAnsTV(x, n, m, lam, max_iter, max_jac_iter, lib_dir="."):
    import numpy.ctypeslib as npct

    array_1d_double = npct.ndpointer(dtype=np.double,
                                     ndim=1, flags='CONTIGUOUS')

    # ipdb.set_trace()
    libl1c = npct.load_library("libl1c", lib_dir)

    libl1c.breg_anistropic_TV.restype = np.int32
    # int breg_anistropic_TV(l1c_int n, l1c_int m, double *uk, double *f,
    #                        double lambda, double tol, l1c_int max_iter,
    #                        l1c_int max_jac_iter)
    libl1c.breg_anistropic_TV.argtypes = [c_int,            # n
                                          c_int,            # m
                                          array_1d_double,  # uk
                                          array_1d_double,  # f
                                          c_double,         # lam
                                          c_double,         # tol
                                          c_int,            # max_iter
                                          c_int]            # max_jac_iter

    uk_out = np.zeros(n*m)
    tol = 0.001
    time0 = time.time()
    stat = libl1c.breg_anistropic_TV(n, m, uk_out, x, lam, tol,
                                     max_iter, max_jac_iter)
    time_total = time.time() - time0
    print("Total python time: %f" % time_total)

    if stat > 0:
        raise Exception("Error in libl1qc.l1qc_dct")

    return uk_out
