clear
clc
fpath = '/home/arnold/matlab/afm-cs/matlab-code/notes/data/cs_sim_CS20NG.mat';
addpath ~/matlab/afm-cs/reconstruction/BP
addpath(genpath('~/matlab/dependencies/SparseLab2.1-Core/Utilities/'))
load(fpath)
%%
whos
tc = matlab.unittest.TestCase.forInteractiveUse;
%%
pix_mask_vec = PixelMatrixToVector(cs_sim.pix_mask);

% y = E*M*x
y_vec = PixelMatrixToVector(cs_sim.Img_sub_sampled);
% y, set of measurements. have to remove all the spots we didn't sample.
y_vec = y_vec(find(pix_mask_vec>0.5));
A = @(x) IDCTfun(x,pix_mask_vec); % E*M
At = @(x) DCTfun(x,pix_mask_vec); %E^T*M^T

% I_bp_dct = l1qc_logbarrier(At(y_vec), A, At, y_vec, 0.1);


x0 = At(y_vec);
b= y_vec;
x = x0;
epsilon = 0.1;
u = (0.95)*abs(x0) + (0.10)*max(abs(x0));
N = length(x0);
% choose initial value of tau so that the duality gap after the first
% step will be about the origial norm
tau = max((2*N+1)/sum(abs(x0)), 1);
cgtol = 1e-8;
cgmaxiter = 200;
cgmaxiter = 200;
lbtol = 1e-3;
newtontol = lbtol;
newtonmaxiter = 50;

[xp, up, ntiter] = l1qc_newton_tmp(x, u, A, At, b, epsilon, tau, newtontol, newtonmaxiter, cgtol, cgmaxiter);
%%
clc
vp = {'AbsTol', 1e-14};

clear l1qc_newton_fcns
[xp2, up2, ntiter2] = l1qc_newton_fcns(x, u, A, At, b, epsilon, tau, newtontol, newtonmaxiter, cgtol, cgmaxiter);

tc.verifyEqual(xp, xp2, vp{:})
tc.verifyEqual(up, up2, vp{:})
%%
dat = loadjson('test_data/line_search_data.json')
r = dat.r;
sum(r)
%%
dat = loadjson('test_data/descent_data.json')