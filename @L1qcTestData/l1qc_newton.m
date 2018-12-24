% l1qc_newton.m
%
% Newton algorithm for log-barrier subproblems for l1 minimization
% with quadratic constraints.
%
% Usage: 
% [xp,up,niter] = l1qc_newton(x0, u0, A, At, b, epsilon, tau, 
%                             newtontol, newtonmaxiter, cgtol, cgmaxiter)
%
% x0,u0 - starting points
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
% tau - Log barrier parameter.
%
% newtontol - Terminate when the Newton decrement is <= newtontol.
%         Default = 1e-3.
%
% newtonmaxiter - Maximum number of iterations.
%         Default = 50.
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


function [xp, up, niter] = l1qc_newton(x0, u0, A, At, b, epsilon, tau, newtontol, newtonmaxiter, cgtol, cgmaxiter, Tii, verbose) 


  % line search parameters
  alpha = 0.01;
  beta = 0.5;
  
  
  % initial point
  x = x0;
  u = u0;
  [r, fu1, fu2, fe, f] =  f_eval(A, x, u, b, epsilon, tau);
  
  
  if verbose
    fprintf('Newton-iter | Functional | Newton decrement |  Stepsize  |  cg-res | cg-iter | backiter |    s     |\n');
  end
  
  for niter=1:newtonmaxiter
    
    % Gradient and Hessian parts.
    atr = At(r);
    ntgz = 1./fu1 - 1./fu2 + 1/fe*atr;
    ntgu = -tau - 1./fu1 - 1./fu2;
    gradf = -(1/tau)*[ntgz; ntgu];
    
    sig11 = 1./fu1.^2 + 1./fu2.^2;
    sig12 = -1./fu1.^2 + 1./fu2.^2;
    sigx = sig11 - sig12.^2./sig11;
    
    w1p = ntgz - sig12./sig11.*ntgu;
    h11pfun = @(z) sigx.*z - (1/fe)*At(A(z)) + 1/fe^2*(atr'*z)*atr;
    [dx, cgres, cgiter] = L1qcTestData.cgsolve(h11pfun, w1p, cgtol, cgmaxiter, 0);
    if (cgres > 1/2)
      disp('Cannot solve system.  Returning previous iterate.  (See Section 4 of notes for more information.)');
      xp = x;  up = u;
      return
    end
    Adx = A(dx);
    du = (1./sig11).*ntgu - (sig12./sig11).*dx;
    
    % minimum step size that stays in the interior
    ifu1 = find((dx-du) > 0); ifu2 = find((-dx-du) > 0);
    aqe = Adx'*Adx;   
    bqe = 2*r'*Adx;   
    cqe = r'*r - epsilon^2;
    root = (-bqe+sqrt(bqe^2-4*aqe*cqe))/(2*aqe);
    if imag(root)~=0
      root = 1;
      fprintf(['Warning: maximum step which satisfies cone constraints has become complex.\n',...
        'Trying smax=1.0 for the cone-constraing portion.\n'])
    end
    smax = min(1,min([...
      -fu1(ifu1)./(dx(ifu1)-du(ifu1));
      -fu2(ifu2)./(-dx(ifu2)-du(ifu2)); ...
      root
      ]));
    s = (0.99)*real(smax);
    
    % backtracking line search
    for backiter=1:32
      xp = x + s*dx;
      up = u + s*du;
      [rp, fu1p, fu2p, fep, fp] =  f_eval(A, xp, up, b, epsilon, tau);
      flin = f + alpha*s*(gradf'*[dx; du]);
      suffdec = (fp <= flin);
      if suffdec
        break
      end
      s = beta*s;
    end
    if ~suffdec
      fprintf(['Stuck on backtracking line search, returning previous iterate.',...
        '(See Section 4 of notes for more information.)\n']);
      xp = x;  
      up = u;
      return
    end
    
    % set up for next iteration
    x = xp; 
    u = up;
    r = rp;
    fu1 = fu1p;  
    fu2 = fu2p;  
    fe = fep; 
    f = fp;
    
    lambda2 = -(gradf'*[dx; du]);
    stepsize = s*norm([dx; du]);
    
    if verbose
      % fprintf('Newton iter |  Functional | Newton decrement | Stepsize | cg-res | backiter|  s \n');
      %            NI         fcnl         dec            sz     cgr       cgI        BI       s
      fprintf('     %3d       %8.3g       %08.3g       % 8.3e   %08.3g     %3d       %2d       %.3g \n',...
        niter, f, lambda2/2, stepsize, cgres, cgiter, backiter, s);
    end
    
    if (lambda2/2 < newtontol)
      break
    end
    
  end
end

function [r, fu1, fu2, fe, f] =  f_eval(A, x, u, b, epsilon, tau)
  % Compute initial point.
  r = A(x) - b;
  fu1 = x - u;
  fu2 = -x - u;
  fe = 1/2*(r'*r - epsilon^2);
  % Should have fe < 0 always. If not, (x,u) is infesible, and so we need the
  % cost to be infinite. But matlab will return a complex number, not Inf or
  % nan. So we force the issue:
  fe = min(fe, 0);
  f = sum(u) - (1/tau)*(sum(log(-fu1)) + sum(log(-fu2)) + log(-fe));

  
end


