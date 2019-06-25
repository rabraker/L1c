#include "config.h"
#include "l1c.h"

/* Helper Routines*/
static int check_pix_idx(l1c_int n, l1c_int *pix_idx, l1c_int max_idx);

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
 * The dimension of the transform `A` is `n` by `(mrow*mcol)`.
 * @param[in]  n Number of elements in pix_idx.
 * @param[in] mrow Number of rows of the underlying signal.
 * @param[in] mcol Number of columns of the underlying signal.
 *             To treat, e.g., an image as a 1D vectorized signal
 *             set `mcol=1` and `mrow = number_of_rows * number_of_columns`.
 * @param[in] dct_mode Whether to setup dct1 (1D) or dct2 (2D)
 *            transforms.
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
int l1c_setup_dct_transforms( l1c_int n, l1c_int mrow, l1c_int mcol, DctMode dct_mode,
                             l1c_int *pix_idx, l1c_AxFuns *ax_funs){

  if (n <=0 || mrow <=0 || mcol <=0){
    return L1C_INVALID_ARGUMENT;
  }

  if (check_pix_idx(n, pix_idx, mrow*mcol - 1)){
    return L1C_INVALID_ARGUMENT;
  }

  if ((mcol == 1) || (mrow == 1) || dct_mode == dct1){
    //call setup_dct1
    return l1c_dct1_setup(n, mrow*mcol, pix_idx, ax_funs);
  }else if ((mcol>1 && mrow > 1) || dct_mode == dct2){
    if (mcol ==1 || mrow == 1){
      return L1C_INCONSISTENT_ARGUMENTS;
    }
    // Call setup_dct2
    return l1c_dct2_setup(n, mrow, mcol, pix_idx, ax_funs);
  }else{
    return L1C_INCONSISTENT_ARGUMENTS;
  }

}

static int check_pix_idx(l1c_int n, l1c_int *pix_idx, l1c_int max_idx) {
  l1c_int idx = 0;
  for (int i = 0; i < n; i++) {
    idx = pix_idx[i];
    if (idx > max_idx || idx < 0) {
      return L1C_INVALID_ARGUMENT;
    }
  }
  return L1C_SUCCESS;
}
