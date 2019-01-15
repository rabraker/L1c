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
% Usage:  xp = l1qc_logbarrier(x0, A, At, b, opts)
%
% x0 - Nx1 vector, initial point.
%
% A - a handle to a function that takes a N vector and returns a K 
%     vector.  
% At - Handle to a function that takes a K vector and returns an N vector.
%
% b - Kx1 vector of observations.
%
%   opts: a struct (create with CsTools.l1qc_opts) with the following fields: 
%         .epsilon : denoising tolerance (as in (1) above). 
%         .mu :      log-barrier weight is multiplied by mu each LB iteration
%                    (try 10).
%         .cgtol:    tolerance for the conjugate gradient (CG) solver (try 1e-8).
%         .cgmaxiter: maximum iterations for the CG solver (try 200).
%         .lbtol :   Determines the total number of log-barrier iterations
%                    (try  1e-3).
%         .newton_tol : newton iterations terminate when the 
%                       newton decrement < newton_tol (try newton_tol = lbtol).
%         .newton_max_iter : Maximum number of newton iterations (try 50).
%         .verbose : How much to print. Silent if verbose =0, prints info each
%         log-barrier iteration if verbose =1, and also each newton iteration
%         if verbose = 2.
%

% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

function xp = l1qc_logbarrier(x0, A, At, b, opts)
  %epsilon, lbtol, mu, cgtol, cgmaxiter, verbose
  if nargin<5
    opts = struct();
  end
  opts = parse_opts(opts);
  
  opts.newtontol = opts.lbtol;
  opts.newtonmaxiter = 50;
  
  N = length(x0);
  
  % starting point --- make sure that it is feasible
  if (norm(A(x0)-b) > opts.epsilon)
    disp('Starting point infeasible; using x0 = At*inv(AAt)*y.');
    AAt = @(z) A(At(z));
    [w, cgres] = cgsolve(AAt, b, opts.cgtol, opts.cgmaxiter, 0);
    if (cgres > 1/2)
      disp('A*At is ill-conditioned: cannot find starting point');
      xp = x0;
      return;
    end
    x0 = At(w);
  end
  x = x0;
  u = (0.95)*abs(x0) + (0.10)*max(abs(x0));
  
  fprintf('Original l1 norm = %.3f, original functional = %.3f\n', sum(abs(x0)), sum(u));
  
  % choose initial value of tau so that the duality gap after the first
  % step will be about the origial norm
  tau = max((2*N+1)/sum(abs(x0)), 1);
  
  lbiter = ceil((log(2*N+1)-log(opts.lbtol)-log(tau))/log(opts.mu));
  fprintf('Number of log barrier iterations = %d\n\n', lbiter);
  
  totaliter = 0;
  total_cg_iter = 0;
  for ii = 1:lbiter
    [xp, up, ntiter, cgiter_ii] = L1qcTestData.l1qc_newton(x, u, A, At, b, opts.epsilon, ...
      tau, opts.newtontol, opts.newtonmaxiter, opts.cgtol, ...
      opts.cgmaxiter, ii, opts.verbose, opts.warm_start_cg);
    totaliter = totaliter + ntiter;
    total_cg_iter = total_cg_iter + cgiter_ii;
    fprintf('\nLog barrier iter = %d, l1 = %.3f, functional = %8.3f, tau = %8.3e, total newton iter = %d, total cg iter: %d\n', ...
      ii, sum(abs(xp)), sum(up), tau, totaliter, total_cg_iter);
    
    x = xp;
    u = up;
    
    tau = opts.mu*tau;
  end


end

function opts = parse_opts(opts)
  flds    = {'lbtol', 'mu', 'cgtol', 'cgmaxiter', 'verbose', 'warm_start_cg'};
  defaults  = [1e-3,    10,    1e-8,    200,        1, 0];
  
  k=1;
  for fld=flds
    if ~isfield(opts, fld{1}) || isempty(opts.(fld{1}))
      opts.(fld{1}) = defaults(k);
    end
    k=k+1;
  end
  
  
end
