function build_cgsolve_test_data(test_data_root)
% Build and save test data for cgsolve.c
% 
% The JSON file will be saved to test_data_root/''cgsolve_small01.json'
%

  cg_test_data_small_path = fullfile(test_data_root, 'cgsolve_small01.json');
  
  A_mat = [4, 1; 1 3];
  b = [1;2];
  
  tol = 1e-6;
  maxiter = 100;
  
  verbose = 1;
  
  cgsolve_tmp(A_mat, b, tol, maxiter, verbose)
  
  
  
  rng(2);
  N = 50;
  A_mat = rand(N,N);
  A_mat = A_mat + A_mat' + 6*eye(N);
  
  min(eig(A_mat))
  
  b = rand(N,1);
  A_fun = @(x) A_mat *x;
  x = cgsolve_local(A_fun, b, tol, maxiter, verbose);
  
  
  A_row = [];
  for i=1:size(A_mat, 1)
    A_row = [A_row, A_mat(i,:)];
  end
  
  data = struct('A', A_row, 'b', b', 'x', x', 'tol', tol, 'max_iter', maxiter);
  
  savejson('', data, cg_test_data_small_path);

end

% cgsolve.m
%
% Solve a symmetric positive definite system Ax = b via conjugate gradients.
%
% Usage: [x, res, iter] = cgsolve(A, b, tol, maxiter, verbose)
%
% A - function handle which computes A_matrix*x.
%
% b - N vector
%
% tol - Desired precision.  Algorithm terminates when 
%    norm(Ax-b)/norm(b) < tol .
%
% maxiter - Maximum number of iterations.
%
% verbose - If 0, do not print out progress messages.
%    If and integer greater than 0, print out progress every 'verbose' iters.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

function [x, res, iter] = cgsolve_local(A, b, tol, maxiter, verbose)

  x = zeros(length(b),1);
  r = b; % = b - A*x0, x0 = 0.
  d = b;
  delta = r'*r;
  delta_init = b'*b;
  
  numiter = 0;
  bestx = x;
  bestres = sqrt(delta/delta_init);
  while ((numiter < maxiter) && (delta > tol^2*delta_init))
    q = A(d);
    alpha = delta/(d'*q);
    x = x + alpha*d;
    
    if (mod(numiter+1,50) == 0)
      r = b - A(x);
    else
      r = r - alpha*q;
    end
    
    deltaold = delta;
    delta = r'*r;
    beta = delta/deltaold;
    d = r + beta*d;
    numiter = numiter + 1;
    if (sqrt(delta/delta_init) < bestres)
      bestx = x;
      bestres = sqrt(delta/delta_init);
    end
    
    if ((verbose) && (mod(numiter,verbose)==0))
      fprintf('cg: Iter = %d, Best residual = %8.3e, Current residual = %8.3e', ...
        numiter, bestres, sqrt(delta/delta_init));
    end
    
  end
  
  if (verbose)
    fprintf('cg: Iterations = %d, best residual = %14.8e', numiter, bestres);
  end
  x = bestx;
  res = bestres;
  iter = numiter;
end