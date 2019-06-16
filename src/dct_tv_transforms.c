#include "config.h"
#include <cblas.h>
#include <stddef.h>
#include <stdio.h>

#include "l1c.h"
#include "TV.h"


#define _JOB_TV_V (1u << 1)
#define _JOB_TV_H (1u << 3)

/* Local functions */
static void destroy();
static void Wv (double *v, double *y);
static void Wtx (double *x, double *v);
static void Rx(double *x, double *y);
static void Rty(double *y, double *x);
static void Av(double *v, double *y);
static void Aty(double *y, double *v);
static void Identity(double *x, double *z);

/* Local variables */
static l1c_int _m;
static l1c_int _mrow;
static l1c_int _mcol;
static l1c_int _p;

static double _alp_h;
static double _alp_v;

static l1c_AxFuns ax_funs_local;
static double *u;
static l1c_int _inc_TVH, _inc_TVV;
static unsigned jobs;

/*
  @TODO This will currently force dct2 if the user wants also TV
  (alp_v>0 or alp_h>0). Fix that, probably have to supply a flag.
  And maybe fix the l1c_setup_dct_transforms as well, for consistency.
*/

int l1c_setup_dctTV_transforms(l1c_int n, l1c_int mrow, l1c_int mcol,
                               double alp_v, double alp_h, DctMode dct_mode,
                               BpMode bp_mode, l1c_int *pix_idx, l1c_AxFuns *ax_funs) {
  jobs = 0;
  int status = L1C_SUCCESS;
  _mrow = mrow;
  _mcol = mcol;
  _m = mrow * mcol;
  _p = _m;
  _alp_v = alp_v;
  _alp_h = alp_h;
  _inc_TVV = 0;
  _inc_TVH = 0;

  ax_funs->norm_W += 1;

  /* Check which jobs we are do to do. We partion local variable u as
     u = [u, u_TVV, u_TVH], for vertical, horizontal parts.
     Also, set inc_TVH, and inc_TVV, to proper increment for aguements below.

  */

  if (alp_v < 0) {
    return L1C_INVALID_ARGUMENT;
  } else if (alp_v > 0) {
    jobs |= _JOB_TV_V;
    _inc_TVV = _m;
  }

  if (alp_h < 0) {
    return L1C_INVALID_ARGUMENT;
  } else if (alp_h > 0) {
    jobs |= _JOB_TV_H;
    _inc_TVH = _inc_TVV + _m;
  }

  if (alp_v > 0 && (mrow <= 2 || mcol <= 2)) {
    return L1C_INVALID_ARGUMENT;
  }

  if (alp_h > 0 && (mrow <= 2 || mcol <= 2)) {
    return L1C_INVALID_ARGUMENT;
  }

  if (alp_v > 0) {
    _p += _m;
    ax_funs->norm_W += 4;
  }
  if (alp_h > 0) {
    ax_funs->norm_W += 4;
    _p += _m;
  }

  ax_funs->n = n;
  ax_funs->m = _m;
  ax_funs->p = _p;

  u = l1c_malloc_double(_p);

  if (!u) {
    return L1C_OUT_OF_MEMORY;
  }

  status = l1c_setup_dct_transforms(n, mrow, mcol, dct_mode,
                                    pix_idx, &ax_funs_local);
  if (status) {
    goto fail;
  }

  ax_funs->Rx = Rx;
  ax_funs->Rty = Rty;

  if (bp_mode == synthesis) {
    ax_funs->Ax = Av;
    ax_funs->Aty = Aty;
    ax_funs->Mx = Wv;
    ax_funs->Mty = Wtx;
    ax_funs->Wz = Identity;
    ax_funs->Wtx = Identity;
  } else {
    ax_funs->Mx = Identity;
    ax_funs->Mty = Identity;
    ax_funs->Ax = Rx;
    ax_funs->Aty = Rty;
    ax_funs->Wz = Wv;
    ax_funs->Wtx = Wtx;
  }
  ax_funs->destroy = destroy;
  ax_funs->data = NULL;

  return status;

 fail:
  l1c_free_double(u);
  return status;
}

static void destroy() {
  // ax_funs_local.destroy();
  l1c_free_double(u);
}

static void Identity(double *x, double *z) {
  cblas_dcopy(_m, x, 1, z, 1);
}

/*
  Implements W * v = [M, Dh, Dv] * [v0;
  v1;
  v2]
  x \in \mathbb{R}^p, p=3*m
*/
static void Wv(double *v, double *y) {

  double *u_TVV = u + _inc_TVV;
  double *u_TVH = u + _inc_TVH;
  ax_funs_local.Mx(v, u);
  cblas_dcopy(_m, u, 1, y, 1);

  if (jobs & _JOB_TV_V) {
    l1c_Dy(_mrow, _mcol, _alp_v, v + _inc_TVV, u_TVV);
    cblas_daxpy(_m, 1, u_TVV, 1, y, 1);
  }

  if (jobs & _JOB_TV_H) {
    l1c_Dx(_mrow, _mcol, _alp_h, v + _inc_TVH, u_TVH);
    cblas_daxpy(_m, 1, u_TVH, 1, y, 1);
  }
}

/*
  Implements W^T * x = [M^T]
  Dh] x
  Dv]

  x \in \mathbb{R}^p, p=3*m
*/
static void Wtx(double *x, double *v) {
  double *v0 = v;
  double *v_dxv = v + _inc_TVV;
  double *v_dxh = v + _inc_TVH;

  ax_funs_local.Mty(x, v0);

  if (jobs & _JOB_TV_V) {
    l1c_DyT(_mrow, _mcol, _alp_v, x, v_dxv);
  }
  if (jobs & _JOB_TV_H) {
    l1c_DxT(_mrow, _mcol, _alp_h, x, v_dxh);
  }
}

static void Rx(double *x, double *y) { ax_funs_local.Rx(x, y); }

static void Rty(double *y, double *x) { ax_funs_local.Rty(y, x); }

static void Av(double *v, double *y) {
  Wv(v, u);
  ax_funs_local.Rx(u, y);
}

static void Aty(double *y, double *v) {

  ax_funs_local.Rty(y, u);
  Wtx(u, v);
}
