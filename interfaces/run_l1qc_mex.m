% This script demonstrates usage of the l1qc interface. We first create an
% image that is a simulation representation of a CS20NG calibration grating,
% (used in Atomic force microscopy). We then create a mu-path mask and sample
% the original image.
%
% Recall the we are solving the problem
%  min_eta ||eta||_1  s.t. ||A*eta - b|| < epsilon.
%
% Here, eta is the estimate of x in the DCT domain, ie, x = M * eta, where M is
% the IDCT transform. The vector b is the sampled data in the standard basis.
%
% This script will likely take a couple minutes to run.
% 
% Unless you are running matlab from the command line, then the l1qc status updates wil 
% will not be instantaneous. Rather, Matlab seems to save up a large buffer
% from mexPrintf and only dumps it to the command window periodically.


clear
addpath ~/matlab/l1c/interfaces

sampling_ratio = 0.1; % take 10% of the pixels
mu_path_len = 40;     % mu-path length in pixels.
N = 512;              % (pixels) Size of square image

% Create the test image.
x_start = 13; % x-offset for the grating holes.
y_stary = 13; 
rng(1);       % Always get the same mask.
X_img_orig = L1qcTestData.make_CS20NG(13,13,N);

% Create the sampling mask. pix_idx is a vector of sampled pixel indeces (which
% assumes the image matrix will be concatenated row-wise). Pix_mask_mat is an
% N by N matrix with 1's where we sampled and zeros elsewhere.

[pix_idx, pix_mask_mat] = L1qcTestData.mu_path_mask(mu_path_len, N, N, sampling_ratio, false);
pix_mask_vec = L1qcTestData.pixmat2vec(pix_mask_mat);

% Put the image matrix into a vector.
img_vec = L1qcTestData.pixmat2vec(X_img_orig);


b = img_vec(pix_idx);      % sub samble the original image.
b = b/max(abs(b)); % normalize to 1

% Define two transformation function handles. If E is the sampling matrix and M
% is the IDCT transform, then A = E * M, and At = M^T * E^T.
At = @(x) L1qcTestData.Atfun_dct(x,pix_idx, length(img_vec)); %E^T*M^T
A = @(x) L1qcTestData.Afun_dct(x,pix_idx); 

% We need the initial guess to l1qc_dct to be feasible, ie, 
% that ||A*eta_0 - b||<epsilon
eta_0 = At(b);



% This will generate an opts struct suitible for passing to l1qc(). 
opts = l1qc_dct_opts('verbose', 2, 'l1_tol', 1e-5);
clear l1qc % In case you just re-compiled.
% system('export MKL_DYNAMIC=TRUE');
% system('export MKL_NUM_THREADS=2');
tic
[eta_est, LBRes]= l1qc_dct(eta_0, b, pix_idx, opts);
tm_mex = toc;
fprintf('mex file time: %.4f\n', tm_mex);

% Put things back in standard basis.
X_est = L1qcTestData.pixvec2mat(idct(eta_est), N);
%return

figure(1); clf
ha = make_axes(gcf);

imagesc(ha(1,1), X_img_orig)
colormap('gray')
title(ha(1,1), 'Original image')

imagesc(ha(1,2), X_img_orig.*pix_mask_mat);
colormap('gray')
title(ha(1,2), '$\mu$-path sampled image')

imagesc(ha(2,1), X_est)
colormap('gray')
mex_title = sprintf('BPDN Reconstruction (matlab)\n %.2f sec, $||\\eta||=%.2f$',...
  tm_mex, norm(eta_est,1));
title(ha(2,1), mex_title)



tic
eta_est_ml = L1qcTestData.l1qc_logbarrier(eta_0, A, At, b, opts);
tm_matlab = toc;

fprintf('matlab l1qc-time: %f\n', tm_matlab);
fprintf('mexl1qc-time:     %f\n', tm_mex);

%%
X_ml = L1qcTestData.pixvec2mat(idct(eta_est_ml), N);
imagesc(ha(2,2), X_ml);
colormap('gray')
ml_title = sprintf('BPDN Reconstruction (matlab)\n %.2f sec, $||\\eta||=%.2f$',...
  tm_matlab, norm(eta_est_ml,1));
title(ha(2,2), ml_title)


function ha = make_axes(fig)
  figure(fig)
  ax1 = axes('Position', [0.1300 0.5482 0.3347 0.3768]);
  ax2 = axes('Position', [0.5703 0.5482 0.3347 0.3768]);
  ax3 = axes('Position', [0.1300 0.0529 0.3347 0.3768]);
  ax4 = axes('Position', [0.5703 0.0457 0.3347 0.3768]);
  ha = [ax1, ax2; ax3, ax4];
end
 