#include "config.h"
#include "l1c_common.h"
#include "dct.h"
#include "dct2.h"


typedef struct XFormPack {
  AxFuns Ax_funs;
  void(*M)(double *x);
  // int(*setup)(l1c_int Nx, l1c_int Ny, l1c_int *pix_mask_idx);
  int(*setup)(l1c_int N, l1c_int M, l1c_int Ny, l1c_int *pix_mask_idx);

  // int(*setup)(l1c_int Nrow, l1c_int Ncol, l1c_int M, l1c_int *pix_idx);
  void(*destroy)(void);

}XFormPack;

void setup_dct_transforms(l1c_int N, l1c_int M, l1c_int Npix, l1c_int *pix_idx, XFormPack *xfms){

  Ax_funs ax_funs;

  if (M = 1){
    //call setup_dct1
    setup_dct(N, Npix, pix_idx);

    ax_funs.Ax = dct_EMx;
    ax_funs.Atx = dct_MtEty;
    ax_funs.AtAx = dct_MtEt_EMx;

    xfms->Ax_funs = ax_funs;
    xfms->destroy = dct_destroy();
  }else if(M>1){
    // Call setup_dct2
    setup_dct2(N, Npix, pix_idx);

    ax_funs.Ax = dct2_EMx;
    ax_funs.Atx = dct2_MtEty;
    ax_funs.AtAx = dct2_MtEt_EMx;

    xfms->Ax_funs = ax_funs;
    xfms->destroy = dct2_destroy();

  }else{
    xfms=NULL;
  }

}

void setup_matrix_transform(l1c_int N, l1c_int M, double *A){

}



/**
   Computes the matrix-vector product y = A * b, for a full matrix A.
   This is a wrapper for cblas_dgemv.
*/
void Ax(l1c_int n, double *x, double *b, void *AX_data){
  double *A = (double *) AX_data;

  cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, A, n, x, 1, 0.0, b, 1);
}
