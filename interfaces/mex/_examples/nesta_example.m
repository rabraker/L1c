% This script creates an image that is a simulation representation of a CS20NG
% calibration grating, (used in Atomic force microscopy). We then create a 
% mu-path mask and sample the original image. We then demonstrate using the
% nesta_dctTV() optimization to reconstruct the image with basis pursuit.

clear
l1c_mex_init_paths();

sampling_ratio = 0.1; % take 10% of the pixels
mu_path_len = 25;     % mu-path length in pixels.
N = 256;              % (pixels) Size of square image

% Create the test image.
try
    rng(1);       % Always get the same mask.
catch ME
    fprintf("CANNOT SET RNG SEED IN OCTAVE\n");
end

X_img_orig = cs20ng_grating(13,13,N);

% Create the sampling mask. pix_idx is a vector of sampled pixel indeces (which
% assumes the image matrix will be concatenated row-wise). Pix_mask_mat is an
% N by N matrix with 1s where we sampled and zeros elsewhere.

[pix_idx, pix_mask_mat] = mu_path_mask(mu_path_len, N, N, sampling_ratio);

% sub samble the original image.
b = X_img_orig(pix_idx);

alpv = .5;
alph = .1;
opts = nesta_opts('alpha_v', alpv,...    % Weight on vertical variation.
                  'alpha_h', alph,...    % Weight on horizontal variation.
                  'verbose', 5,...       % print every 5th iteration.
                  'dct_mode', 1, ...     % use 1d dct
                  'bp_mode', 'analysis',... %Can also select synthesis mode.
                  'tol', 1e-5,...        % Read the paper...
                  'mu', 1e-5);           % Read the paper...

[x_est, status] = nesta_dctTV(N, N, b, pix_idx, opts);

figure(1, 'Position', [0, 0 1400, 450])
subplot(1,3,1)
imagesc(X_img_orig)
colormap('gray')
title('Original')

subplot(1,3,2)
imagesc(pix_mask_mat)
colormap('gray')
title('Original')

subplot(1,3,3)
imagesc(x_est)
colormap('gray')
title('Reconstruction')
