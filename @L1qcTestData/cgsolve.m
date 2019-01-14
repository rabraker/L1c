function [x, res, iter] = cgsolve(A, b, tol, maxiter, verbose, x0)
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
% x0 - initial guess.
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

  if nargin < 6 || isempty(x0)
    x = zeros(length(b), 1);
  else
    x = x0;
  end
  

  r = b - A(x);
  p = r;
  
  delta = r'*r;
  delta_init = b'*b;
  
  bestx = x;
  bestres = sqrt(delta/delta_init);
  if verbose
    fprintf('cg: |Iter| Best resid | Current resid| alpha | beta   |   delta  |\n');
  end
 for iter = 1:maxiter
    q = A(p);
    alpha = delta/(p'*q);
    x = x + alpha*p;

%     if (mod(iter+1,50) == 0)
%       r = b - A(x);
%       p = r;
%       delta = r'*r;
%       continue;
%     end

    r = r - alpha*q;

    deltaold = delta;
    delta = r'*r;
    beta = delta/deltaold;
    p = r + beta*p;
    
    rel_res = sqrt(delta/delta_init); % sqrt for norm.
    if rel_res < bestres
      bestx = x;
      bestres = rel_res;
    end
    
    if ((verbose) && (mod(iter,verbose)==0))
      %         iter              br             cr    alpha beta delta
      fprintf('  %d,   %.16e, %.16e, %.16e, %.16e, %.16e  \n', ...
        iter, bestres, rel_res, alpha, beta, delta);
    end
    if (rel_res < tol)
      break;
    end
    
  end
  
  if (verbose)
    fprintf('cg: Iterations = %d, best residual = %14.8e\n', iter, bestres);
  end
  x = bestx;
  res = bestres;
  
end