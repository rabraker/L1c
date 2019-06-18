#!/usr/bin/env python3

import numpy as np
import L1cTestDataUtils as TDU


def make_CS20NG(npix, x_start=10, y_start=10, nholes=10):

    pitch = npix/nholes
    hole_width = npix/2/nholes

    m_pitch = np.int32(np.ceil(pitch))
    n_hole = np.int32(np.ceil(hole_width))
    img_mat = np.ones((npix, npix))

    # First lets, create a prototype hole. Then we will insert it into
    # the master image periodically.

    sq = np.ones((n_hole, n_hole))

    # circle radius
    r = n_hole/2

    # center
    xc = n_hole/2
    yc = n_hole/2

    for x_pt in range(0, n_hole):
        for y_pt in range(0, n_hole):
            rad_pt = np.sqrt((x_pt - xc)**2 + (y_pt - yc)**2)

            if rad_pt <= r:
                sq[y_pt, x_pt] = 0  # make it black

    x_lft = x_start
    n, m = sq.shape
    while x_lft+n_hole < npix + n_hole:
        y_tp = y_start
        x_stride = max(min(sq.shape[1], npix-x_lft), 0)
        while y_tp+n_hole <= npix + n_hole:
            y_stride = max(min(n, npix-y_tp), 0)
            img_mat[y_tp:y_tp+y_stride, x_lft:x_lft+x_stride] = sq[0:y_stride, 0:x_stride]

            y_tp = y_tp + m_pitch
        x_lft = x_lft + m_pitch
        # img_mat = img_mat(1:npix, 1:npix)
    return img_mat


def mu_path_mask(mupath_len, npix, samplingRatio, repeat_sampling=False):
    """
 Return Mu path sampling mask

 Arguments
    ----------
    mupath_len : size (in pixels) of mu path pattern
    npix : image size
    samplingRatio : (Percentage) - total pixels to sample = samplingRatio*n*m
    repeat_sampling : (true|false) Not implemented. If true, repeat sampling
                      allowed.If false: repeat sampling not allowed.

    Returns
    --------
    pix_idx : a vector of sampled pixel indeces.
    pix_mask_mat : an n by m mu path pattern mask, which contains 1's in the
                  mu-path areas, and zeros elsewhere.
    """

    pix_mask_mat = np.zeros((npix, npix))

    while np.sum(pix_mask_mat) < samplingRatio * npix**2:
        rand_i = np.random.randint(npix)
        rand_j = np.random.randint(2-mupath_len, high=npix)
        x_start = max(1, rand_j)
        x_end = min(rand_j+mupath_len-1, npix)
        pix_mask_mat[rand_i, x_start:x_end] = 1

    pix_idx, = np.where(pix_mask_mat.flatten() > 0.5)
    return pix_idx, pix_mask_mat


def build_cs20ng_test_data(npix, data_path=None):
    img = make_CS20NG(npix)
    pix_idx, pix_mask_mat = mu_path_mask(20, npix, 0.1)

    data = {'x_orig': img,
            'mtot': img.flatten().shape[0],
            'mrow': img.shape[0],
            'mcol': img.shape[1],
            'b': img.flatten()[pix_idx],
            'pix_idx': pix_idx,
            'one_based_index': 0}

    data = TDU.jsonify(data)

    if data_path is not None:
        TDU.save_json(data, data_path)

    return data


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        test_data_path = TDU.data_dir()
        image_data_path_256 = test_data_path+"/example_img_data.json"
        build_cs20ng_test_data(256, data_path=image_data_path_256)

        # Regression Test: This size caused segfaults, due to alignment.
        image_data_path_127 = test_data_path+"/example_img_data127.json"
        build_cs20ng_test_data(127, data_path=image_data_path_127)
    else:
        image_data_path = sys.argv[1]
        build_cs20ng_test_data(256, data_path=image_data_path)



# plt.figure(1)
# plt.imshow(img)

# plt.figure(2)
# plt.imshow(pix_mask_mat)

# plt.show()
