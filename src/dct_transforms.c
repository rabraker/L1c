#include "config.h"
#include "l1c.h"

/** @file
   The high level interface to the DCT based transforms, suitible
   for passing to, e.g., l1qc_newton().
 */


/**
 * @ingroup transforms
 *
 * Purpose
 * =======
 * This is the high-level interface to dct1_setup() and dct2_setup()
 * The function will initialize the data structures for the sub-sampled
 * DCT1/DCT2 transformations. If `mx > 1`, dct2_setup() will called. If
 * `mx == 1`, dct1_setup() will be called. Otherwise, the function returns
 * `L1C_DCT_INIT_FAILURE`.
 *
 *
 *
 * @param[in] Nx Number of rows of the underlying signal.
 * @param[in] Mx Number of columns of the underlying signal.
 *             To treat, e.g., an image as a 1D vectorized signal
 *             set Mx=1 and Nx = number_of_rows * number_of_columns.
 * @param[in]  Ny Number of elements in pix_idx.
 * @param[in]  pix_idx indeces of locations of the subsampling. For both
 *             DCT1 and DCT2, this vector should be the same.
 * @param[out] ax_funs A structure of function pointers which will be populated.
 *             On successfull exit, The fields Ax, Aty, AtAx, destroy, and M
 *             will be non-null. The fields MT, E, and ET will be null.
 *
 *
 * @return      0 if succesfull. If unsuccesfull, returns L1C_DCT_INIT_FAILURE
 *              or L1C_OUT_OF_MEMORY.
 *
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
