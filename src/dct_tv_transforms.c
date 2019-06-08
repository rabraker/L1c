#include "config.h"
#include "l1c.h"
#include "TV.h"
#include <stddef.h>


/* Local functions */
static void destroy();
static void Wv (double *v, double *y);
static void Wtx (double *x, double *v);
static void Ex(double *x, double *y);
static void Ety(double *y, double *x);
static void Av(double *v, double *y);
static void Aty(double *y, double *v);


/* Local variables */
static l1c_int _m;
static l1c_int _mrow;
static l1c_int _mcol;
static l1c_int _p;

static double _alp_h;
static double _alp_v;

static l1c_AxFuns ax_funs_local;
static double *u;




/*
  @TODO This will currently force dct2 if the user wants also TV
  (alp_v>0 or alp_h>0). Fix that, probably have to supply a flag.
  And maybe fix the l1c_setup_dct_transforms as well, for consistency.
 */

int l1c_setup_dctTV_transforms(l1c_int n, l1c_int mrow, l1c_int mcol,
                               double alp_v, double alp_h, l1c_int *pix_idx,
                               l1c_AxFuns *ax_funs){
  int status = 0;
  _mrow = mrow;
  _mcol = mcol;
  _m = mrow * mcol;

  _p = _m;
  ax_funs->norm_M += 1;

  if (alp_v <0){
    return L1C_INVALID_ARGUMENT;
  }
  if (alp_h <0){
    return L1C_INVALID_ARGUMENT;
  }

  if (alp_v > 0 && (mrow <=2 || mcol <=2) ){
    return L1C_INVALID_ARGUMENT;
  }
  if (alp_h > 0 && (mrow <=2 || mcol <=2) ){
    return L1C_INVALID_ARGUMENT;
  }


  if (alp_v > 0){
    _p += _m;
    ax_funs->norm_M += 4;
  }
  if (alp_h > 0) {
    ax_funs->norm_M += 4;
    _p += _m;
  }

  u = l1c_malloc_double(_p);

  if (!u){
    return L1C_OUT_OF_MEMORY;
  }

  status =l1c_setup_dct_transforms(mrow, mcol, n, pix_idx, &ax_funs_local);
  if (status){
    goto fail;
  }

  ax_funs->Mx = Wv;
  ax_funs->Mty = Wtx;
  ax_funs->Ex = Ex;
  ax_funs->Ety = Ety;
  ax_funs->Ax = Av;
  ax_funs->Aty = Aty;
  ax_funs->destroy = destroy;
  ax_funs->data = NULL;

  return status;

 fail:
  l1c_free_double(u);
  return status;
}


static void destroy(){
  ax_funs_local.destroy();
  l1c_free_double(u);
}

/*
  Implements W * v = [M, Dh, Dv] * [v0;
                                    v1;
                                    v2]
  x \in \mathbb{R}^p, p=3*m
 */
static void Wv (double *v, double *y){
  double *u0 = u;
  double *u1 = u + _m;
  double *u2 = u + 2*_m;

  ax_funs_local.Mx(v, u0);
  l1c_Dx(_mrow, _mcol, _alp_h, v+_m, u1);
  l1c_Dy(_mrow, _mcol, _alp_v, v+2*_m, u2);
  // l1c_Dx(_mrow, _mcol, v+_m, u1);
  // l1c_Dy(_mrow, _mcol, v+2*_m, u2);

  for (int i=0; i<_m; i++){
    y[i] = u0[i] + u1[i] + u[2];
  }
}

/*
  Implements W^T * x = [M^T]
                        Dh] x
                        Dv]

 x \in \mathbb{R}^p, p=3*m
*/
static void Wtx (double *x, double *v){

  ax_funs_local.Mty(x, v);
  l1c_DxT(_mrow, _mcol, _alp_h, x, v + _m);
  l1c_DyT(_mrow, _mcol, _alp_v, x, v + 2*_m);

}

static void Ex(double *x, double *y){
  ax_funs_local.Ex(x, y);
}

static void Ety(double *y, double *x){
  ax_funs_local.Ety(y, x);
}

static void Av(double *v, double *y){
  Wv(v, u);
  ax_funs_local.Ex(u, y);
}

static void Aty(double *y, double *v){

  ax_funs_local.Ety(y, u);
  Wtx(u, v);
}
