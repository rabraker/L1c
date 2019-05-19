% This script creates an image that is a simulation representation of a CS20NG
% calibration grating, (used in Atomic force microscopy). We then create a 
% mu-path mask and sample the original image. We then demonstrate using the
% l1qc_newton() optimization to reconstruct the image with basis pursuit.
% Finally, we demonstrate the anistrpic TV denoising, to remove some noise from
% the BP reconstruction.

clear

addpath ../../interfaces/mex
% This needs to point to the current build directory, ie whereever the mex files
% ended up.
addpath ../../build/interfaces/mex

sampling_ratio = 0.1; % take 10% of the pixels
mu_path_len = 25;     % mu-path length in pixels.
N = 256;              % (pixels) Size of square image

% Create the test image.
x_start = 13; % x-offset for the grating holes.
y_stary = 13; 
rng(1);       % Always get the same mask.
X_img_orig = cs20ng_grating(13,13,N);

% Create the sampling mask. pix_idx is a vector of sampled pixel indeces (which
% assumes the image matrix will be concatenated row-wise). Pix_mask_mat is an
% N by N matrix with 1s where we sampled and zeros elsewhere.

[pix_idx, pix_mask_mat] = mu_path_mask(mu_path_len, N, N, sampling_ratio);

% Put the image matrix into a vector.
img_vec = reshape(X_img_orig', N*N,1);

b = img_vec(pix_idx);      % sub samble the original image.
b = b/max(abs(b)); % normalize to 1

n = length(img_vec)
opts = l1qc_dct_opts('verbose', 2, 'l1_tol', 1e-5);

[x_est, LBRes]= l1qc_dct(n, 1, b, pix_idx, opts);

X_bp = reshape(x_est, N, N)';

figure(1)
subplot(2, 2, 1)
imagesc(X_img_orig)
colormap('gray')

subplot(2, 2, 2)
imagesc(pix_mask_mat)
colormap('gray')

subplot(2, 2, 3)
imagesc(X_bp)
colormap('gray')

% We can help with noise in the BP reconstruction with TV denoising.
mu = 5;
max_iter = 1000;
tol = 0.001;
X_tv = breg_anistropic_TV(X_bp, mu, tol, max_iter);

subplot(2, 2, 4)
imagesc(X_tv)
colormap('gray')



