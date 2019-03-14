% Uopt = breg_anistropic_TV_mex(Img, mu, tol, max_iter)
%
%  Solve the anistropic TV denoising problem
% 
%  min ||\\nabla_x u||_1 + ||\\nabla_y||_1 + 0.5\\mu ||u - Img||_2 \n      (1)
%   u
% 
%   using Bregman Splitting. This algorithm is developed in the paper 
%   "The Split Bregman Method for L1-Regularized Problems," 
%   by Tom Goldstein and Stanley Osher. 
% 
%  Arguments
%  ------
%   Img : (double, m by n) Noisy input matrix. Must have two dimensions and be
%          at least three by three.
%   mu :  (double scalar) see (1). Optional, default 5.
%   tol:  (double scalar) Optional, default 0.001. The optimization stops when
%          ||u_k - u_{k-1}||/||u_{k-1}|| < tol
% 
%   max_iter : (int scalar) max number of iterations. Optional, default 1000.
% 
% 
%  Returns
%  -------
%    Uopt : (double m by n matrix).
% 
% 