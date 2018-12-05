% Create a small set of test data
N = 50;

Ts = 1/50;
To = 0.25;

omega =2*pi/To;

k = [ 0:N-1]';

x = sin(omega * Ts *k);

pix_idx = [2, 10, 15, 20, 25, 30, 35, 40, 45, 49];

Ny = length(pix_idx);
pix_mask = zeros(N,1);
pix_mask(pix_idx) = 1;

yy = IDCTfun(x, pix_mask);

% account for scaling difference
% f_AtA = f_AtA * 2*sqrt(N/2)
% f_AtA(1) = f_AtA(1) * sqrt(2);

data = struct('x0', x(:)', 'x1', yy(:)', 'pix_idx', pix_idx(:)'-1);
savejson('', data, 'test_data/dct_small_EMx.json');

%
xx = DCTfun(yy, pix_mask);
xx_save = xx; 
xx_save(1) = xx_save(1) * sqrt(2); % matlabs weird normalization.

data = struct('x0', xx_save(:)', 'x1', yy(:)', 'pix_idx', pix_idx(:)'-1);
savejson('', data, 'test_data/dct_small_MtEty.json');
%
data = struct('x0', x(:)', 'x1', xx_save(:)', 'pix_idx', pix_idx(:)'-1);
savejson('', data, 'test_data/dct_small_MtEt_EMx.json');

%%
% N = 50;
N = 512*512;

Ts = 1/50
To = 0.25

omega =2*pi/To

k = [ 0:N-1]';

x = sin(omega * Ts *k);

% plot(Ts*k, y)
clc
[x, dct(x) * 2 *sqrt(50/2)];

pix_idx = [1, 10, 15, 20, 25, 30, 35, 40, 45, 49]+1;
% pix_idx = 1:50;
pix_idx = 1:23617;

% Ny = length(pix_idx);
% pix_mask(pix_idx) = 1;
% pix_mask = zeros(N,1);
pix_mask = zeros(N,1);
pix_mask(pix_idx) = 1;

clc

yy = IDCTfun(x, pix_mask);
f_AtA = DCTfun(yy, pix_mask);

% account for scaling difference
% f_AtA = f_AtA * 2*sqrt(N/2)
% f_AtA(1) = f_AtA(1) * sqrt(2);

% [f_AtA, y]


clc

y = IDCTfun(x, pix_mask);
x = DCTfun(y, pix_mask);

x(1: min(25, length(y)))
%%
tic

%%
x = fun_AtA(x, pix_mask);
x(1:25)
toc/1


function x = fun_AtA(x, pix_mask)
  
  for k=1:1
    y = IDCTfun(x, pix_mask);
    x = DCTfun(y, pix_mask);
  end
  
end