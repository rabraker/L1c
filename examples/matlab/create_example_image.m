% This script creates an image that is a simulation representation of a CS20NG
% calibration grating, (used in Atomic force microscopy). We then create a 
% mu-path mask and sample the original image. The data is saved as a json file,
% so it can be loaded by the various example/test scripts.

clear
addpath ~/matlab/l1c/interfaces
addpath ~/matlab/l1c/

sampling_ratio = 0.1; % take 10% of the pixels
mu_path_len = 40;     % mu-path length in pixels.
N = 256;              % (pixels) Size of square image

% Create the test image.
x_start = 13; % x-offset for the grating holes.
y_stary = 13; 
rng(1);       % Always get the same mask.
X_img_orig = L1qcTestData.make_CS20NG(13,13,N);

% Create the sampling mask. pix_idx is a vector of sampled pixel indeces (which
% assumes the image matrix will be concatenated row-wise). Pix_mask_mat is an
% N by N matrix with 1s where we sampled and zeros elsewhere.

[pix_idx, pix_mask_mat] = L1qcTestData.mu_path_mask(mu_path_len, N, N, sampling_ratio, false);
pix_mask_vec = L1qcTestData.pixmat2vec(pix_mask_mat);

% Put the image matrix into a vector.
img_vec = L1qcTestData.pixmat2vec(X_img_orig);

b = img_vec(pix_idx);      % sub samble the original image.
b = b/max(abs(b)); % normalize to 1

n = length(img_vec)
opts = l1qc_dct_opts('verbose', 2, 'l1_tol', 1e-5);

[x_est, lb_res] = l1qc_so_wrap(n, b, pix_idx, opts);


opts.FileName = 'example_img_data.json';
opts.FloatFormat = '%.15f';
savejson('', struct('x_orig', img_vec(:)', 'N', length(img_vec),...
                    'b', b(:)', 'pix_idx', pix_idx(:)'-1), opts);
