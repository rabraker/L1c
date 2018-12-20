% l1qc_logbarrier.m
%
% Solve quadratically constrained l1 minimization:
% min ||x||_1   s.t.  ||Ax - b||_2 <= \epsilon
%
% Reformulate as the second-order cone program
% min_{x,u}  sum(u)   s.t.    x - u <= 0,
%                            -x - u <= 0,
%      1/2(||Ax-b||^2 - \epsilon^2) <= 0
% and use a log barrier algorithm.
%
% Usage:  xp = l1qc_logbarrier(x0, A, At, b, epsilon, lbtol, mu, cgtol, cgmaxiter)
%
% x0 - Nx1 vector, initial point.
%
% A - Either a handle to a function that takes a N vector and returns a K 
%     vector , or a KxN matrix.  If A is a function handle, the algorithm
%     operates in "largescale" mode, solving the Newton systems via the
%     Conjugate Gradients algorithm.
%
% At - Handle to a function that takes a K vector and returns an N vector.
%      If A is a KxN matrix, At is ignored.
%
% b - Kx1 vector of observations.
%
% epsilon - scalar, constraint relaxation parameter
%
% lbtol - The log barrier algorithm terminates when the duality gap <= lbtol.
%         Also, the number of log barrier iterations is completely
%         determined by lbtol.
%         Default = 1e-3.
%
% mu - Factor by which to increase the barrier constant at each iteration.
%      Default = 10.
%
% cgtol - Tolerance for Conjugate Gradients; ignored if A is a matrix.
%     Default = 1e-8.
%
% cgmaxiter - Maximum number of iterations for Conjugate Gradients; ignored
%     if A is a matrix.
%     Default = 200.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

function xp = build_logbarrier_test_data(data_root, lbiter)  
  fpath = '/home/arnold/matlab/afm-cs/matlab-code/notes/data/cs_sim_CS20NG.mat';
  addpath ~/matlab/afm-cs/matlab-code/functions
  addpath ~/matlab/afm-cs/reconstruction/BP
  addpath ~/matlab/afm-cs/reconstruction/BP/l1magic/Optimization/
  
  dat = load(fpath);
  cs_sim = dat.cs_sim;
  pix_mask_vec = PixelMatrixToVector(cs_sim.pix_mask);
  pix_idx = find(pix_mask_vec>0.5);
  
  b = PixelMatrixToVector(cs_sim.Img_original);
  b = b/max(b); % normalize to 1.
  b = b(pix_idx); 
  
  
  A = @(x) IDCTfun(x,pix_mask_vec);
  At = @(x) DCTfun(x,pix_mask_vec);
  x0 = At(b);
 
  epsilon = 0.1;
  lbtol = 1e-3;
  mu = 10;
  cgtol = 1e-9;
  cgmaxiter = 400;

  newtontol = lbtol;
  newtonmaxiter = 50;

  N = length(x0);

  % starting point --- make sure that it is feasible
  if (norm(A(x0)-b) > epsilon)
    error('Starting point infeasible; using x0 = At*inv(AAt)*y.');
%     AAt = @(z) A(At(z));
%     [w, cgres] = cgsolve(AAt, b, cgtol, cgmaxiter, 0);
%     if (cgres > 1/2)
%       disp('A*At is ill-conditioned: cannot find starting point');
%       xp = x0;
%       return;
%     end
%     x0 = At(w);
  end
  x = x0;
  u = (0.95)*abs(x0) + (0.10)*max(abs(x0));

  fprintf('Original l1 norm = %.3f, original functional = %.3f\n', sum(abs(x0)), sum(u));

  % choose initial value of tau so that the duality gap after the first
  % step will be about the origial norm
  tau = max((2*N+1)/sum(abs(x0)), 1);
  
  if nargin <2
    lbiter = ceil((log(2*N+1)-log(lbtol)-log(tau))/log(mu));
  end
  
  fprintf('Number of log barrier iterations = %d\n\n', lbiter);

  totaliter = 0;
  jopts.FloatFormat = '%.20f';
  for ii = 1:lbiter
   
    [xp, up, ntiter] = L1qcTestData.l1qc_newton(x, u, A, At, b, epsilon, tau,...
      newtontol, newtonmaxiter, cgtol, cgmaxiter, ii);
    totaliter = totaliter + ntiter;
    
    fprintf(['\nLog barrier iter = %d, l1 = %.3f, functional = %8.3f, ',...
      'tau = %8.3e, total newton iter = %d\n'], ...
            ii, sum(abs(xp)), sum(up), tau, totaliter);
    
    fname_iter = sprintf('lb_test_data_iter_%d.json', ii);
    jopts.FileName = fullfile(data_root, fname_iter);
    savejson('', struct('x0', x(:)', 'b', b(:)',...
      'epsilon', epsilon, 'tau', tau, 'lbiter', lbiter, 'lbtol', lbtol,...
      'newtontol', newtontol, 'newtonmaxiter', newtonmaxiter,...
      'cgtol', cgtol, 'cgmaxiter', cgmaxiter, 'mu', mu, 'pix_idx', pix_idx(:)'-1,...
      'xp', xp(:)', 'up', up(:)', 'tau_next', mu*tau), jopts);
    
    tau = mu*tau;
    x = xp;
    u = up;
    
  end

end