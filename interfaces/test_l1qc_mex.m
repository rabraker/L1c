fpath = '/home/arnold/matlab/afm-cs/matlab-code/notes/data/cs_sim_CS20NG.mat';
addpath ~/matlab/afm-cs/matlab-code/functions
addpath ~/matlab/afm-cs/reconstruction/BP/

dat = load(fpath);
cs_sim = dat.cs_sim;


pix_mask_vec = PixelMatrixToVector(cs_sim.pix_mask);

pix_idx = find(pix_mask_vec>0.5);
img_vec = L1qcTestData.pixmat2vec(cs_sim.Img_sub_sampled);
% y = E*M*x

% y, set of measurements. have to remove all the spots we didn't sample.
y_vec = img_vec(pix_idx);
y_vec = y_vec/max(y_vec); % normalize to 1
At = @(x) L1qcTestData.Atfun_dct(x,pix_idx, length(img_vec)); %E^T*M^T



x0 = At(y_vec);
b= y_vec;
x = x0;

clear l1qc
opts.epsilon = 0.1;
opts.mu = 10;
opts.cgtol = 1e-8;
opts.cgmaxiter = 200;
opts.lbtol = 1e-3;
opts.newton_tol = opts.lbtol;
opts.newton_max_iter = 50;
opts.verbose = 2;
opts.warm_start_cg = 0;
% system('export MKL_DYNAMIC=TRUE');
% system('export MKL_NUM_THREADS=2');
tic
[ximg, LBRes]= l1qc(x0, b, pix_idx-1, opts);
tm_c = toc;
fprintf('mex file time: %.4f\n', tm_c);

X = PixelVectorToMatrix(idct(ximg), [512, 512]);
%%
figure(1)
subplot(2,2, [1,3])
imagesc(X);
colormap('gray')
%%
x0 = At(y_vec);
tic
ximg_ml = L1qcTestData.l1qc_logbarrier(x0, A, At, y_vec, opts.epsilon, opts.lbtol, opts.mu,...
  opts.cgtol, opts.cgmaxiter, opts.verbose);
tml = toc;

fprintf('matlab l1qc-time: %f\n', tml);
fprintf('mexl1qc-time:     %f\n', tm_c);

subplot(2,2, [2,4])
X_ml = PixelVectorToMatrix(idct(ximg_ml), [512, 512]);
imagesc(X_ml);
colormap('gray')
 