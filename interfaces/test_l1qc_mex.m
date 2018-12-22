fpath = '/home/arnold/matlab/afm-cs/matlab-code/notes/data/cs_sim_CS20NG.mat';
addpath ~/matlab/afm-cs/matlab-code/functions
addpath ~/matlab/afm-cs/reconstruction/BP/

dat = load(fpath);
cs_sim = dat.cs_sim;


pix_mask_vec = PixelMatrixToVector(cs_sim.pix_mask);

pix_idx = find(pix_mask_vec>0.5);

% y = E*M*x
y_vec = PixelMatrixToVector(cs_sim.Img_sub_sampled);
% y, set of measurements. have to remove all the spots we didn't sample.
y_vec = y_vec(pix_mask_vec>0.5);
y_vec = y_vec/max(y_vec); % normalize to 1
A = @(x) L1qcTestData.IDCTfun(x,pix_mask_vec); % E*M
At = @(x) L1qcTestData.DCTfun(x,pix_mask_vec); %E^T*M^T



x0 = At(y_vec);
b= y_vec;
x = x0;
opts.epsilon = 0.1;
u = (0.95)*abs(x0) + (0.10)*max(abs(x0));
N = length(x0);
% choose initial value of tau so that the duality gap after the first
% step will be about the origial norm

clear l1qc
opts.mu = 10;
opts.cgtol = 1e-8;
opts.cgmaxiter = 200;
opts.lbtol = 1e-3;
opts.newton_tol = opts.lbtol;
opts.newton_max_iter = 50;
opts.verbose = 1;

tic
ximg= l1qc(x0, b, pix_idx-1, opts);
tm_c = toc;

X = PixelVectorToMatrix(idct(ximg), [512, 512]);

figure(1)
subplot(2,2, [1,3])
imagesc(X);
colormap('gray')
 %%
 tic
 ximg_ml = l1qc_logbarrier_tmp(x0, A, At, b, opts.epsilon, opts.lbtol, opts.mu,...
   opts.cgtol, opts.cgmaxiter, opts.verbose);
 tml = toc;
 %%
 fprintf('matlab l1qc-time: %f\n', tml);
  fprintf('mexl1qc-time:     %f\n', tm_c);
  
 subplot(2,2, [2,4])
 X_ml = PixelVectorToMatrix(idct(ximg_ml), [512, 512]);
 
 