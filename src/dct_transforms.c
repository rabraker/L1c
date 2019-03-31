#include "config.h"
#include "l1c.h"

/*
  This is the high-level interface to dct1.c and dct2.c
 */



int l1c_setup_dct_transforms(l1c_int Nx, l1c_int Mx, l1c_int Ny,
                             l1c_int *pix_idx, l1c_AxFuns *ax_funs){

  if (Mx == 1){
    //call setup_dct1
    return l1c_dct1_setup(Nx, Ny, pix_idx, ax_funs);
  }else if(Mx>1){
    // Call setup_dct2
    return l1c_dct2_setup(Nx, Mx, Ny, pix_idx, ax_funs);

  }else{
    return L1C_DCT_INIT_FAILURE;
  }

}
