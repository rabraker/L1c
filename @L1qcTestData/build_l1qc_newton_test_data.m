function build_l1qc_newton_test_data(test_data_root)
clear
clc
fpath = '/home/arnold/matlab/afm-cs/matlab-code/notes/data/cs_sim_CS20NG.mat';
addpath ~/matlab/afm-cs/reconstruction/BP
addpath(genpath('~/matlab/dependencies/SparseLab2.1-Core/Utilities/'))
dat = load(fpath);
cs_sim = dat.cs_sim;


pix_mask_vec = PixelMatrixToVector(cs_sim.pix_mask);

% y = E*M*x
y_vec = PixelMatrixToVector(cs_sim.Img_sub_sampled);
% y, set of measurements. have to remove all the spots we didn't sample.
y_vec = y_vec(find(pix_mask_vec>0.5));
A = @(x) IDCTfun(x,pix_mask_vec); % E*M
At = @(x) DCTfun(x,pix_mask_vec); %E^T*M^T

x0 = At(y_vec);
b= y_vec;
x = x0;
opts.epsilon = 0.1;
u = (0.95)*abs(x0) + (0.10)*max(abs(x0));
N = length(x0);
% choose initial value of tau so that the duality gap after the first
% step will be about the origial norm
opts.tau = max((2*N+1)/sum(abs(x0)), 1);
opts.cgtol = 1e-8;
opts.cgmaxiter = 200;
opts.lbtol = 1e-3;
opts.newtontol = lbtol;
opts.newtonmaxiter = 50;


vp = {'AbsTol', 1e-14};


[xp2, up2, ntiter2] = l1qc_newton_local(x, u, A, At, b, opts, test_data_root);
  
end


function [x, u, niter] = l1qc_newton_fcns(x0, u0, A, At, b, opts, data_root) 
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


epsilon = opts.epsilon;
tau = opts.tau;
newtontol = opts.newtontol;
newtonmaxiter = opts.newtonmaxiter;
cgtol = opts.cgtol;
cgmaxiter = opts.cgmaxiter;

  % line search parameters
  alpha = 0.01;
  beta = 0.5;  

  % initial point
  x = x0;
  u = u0;
  r = A(x) - b; 

  fu1 = x - u;
  fu2 = -x - u;
  fe = 1/2*(r'*r - epsilon^2);
  f = sum(u) - (1/tau)*(sum(log(-fu1)) + sum(log(-fu2)) + log(-fe));

  jopts.FloatFormat ='%.15g';
  for niter=1:newtonmaxiter
    [dx, du, gradf, cgres, gd] = compute_descent(fu1, fu2, r, fe, tau, ...
                                             cgtol, cgmaxiter, A, At);

    if niter == 5
      jopts.FileName = fullfile(data_root, 'descent_data.json');
    savejson('', struct('dx', dx(:)', 'du', du(:)', 'gradf', gradf(:)', 'cgres', cgres, 'fu1', ...
                        fu1(:)', 'fu2', fu2(:)', 'r', r(:)', 'atr', gd.atr(:)',...
                        'ntgu', gd.ntgu(:)', 'sig11', gd.sig11(:)', 'sig12', gd.sig12(:)',...
                        'w1p', gd.w1p(:)', 'sigx', gd.sigx(:)', 'fe', fe, 'tau', ...
                        tau, 'cgtol', cgtol, 'cgmaxiter', ...
                        cgmaxiter), jopts);
    end
   
    if (cgres > 1/2)
      disp('Cannot solve system.  Returning previous iterate.  (See Section 4 of notes for more information.)');
      xp = x;  up = u;
      return
    end
    
    
    Adx = A(dx);
    
    % get the maximum step size.
    s = find_max_step(dx, du, Adx, fu1, fu2, r, epsilon);
    if niter == 5
      jopts.FileName = fullfile(data_root,'max_stepsize_data.json');
      savejson('', struct('s', s, 'dx', dx(:)', 'Adx', Adx(:)', 'fu1', fu1(:)', ...
                          'fu2', fu2(:)', 'r', r(:)', 'epsilon', epsilon), jopts);
    end
    
    
    [xp, up, rp, fu1p, fu2p, fep, fp, BI] = line_search(x, u, r, dx, du, Adx, gradf, f, tau, ...
                                    epsilon, alpha, beta, s);
    if niter == 5
      jopts.FileName = fullfile(data_root, 'line_search_data.json');
      savejson('', struct('xp', xp(:)', 'up', up(:)', 'rp', rp(:)', 'fu1p', fu1p(:)', ...
                          'fu2p', fu2p(:)', 'fep', fep(:)', 'fp', fp(:)', 'x', ...
                          x(:)', 'u', u(:)', 'r', r(:)', 'dx', dx(:)', 'du', du(:)', ...
                          'Adx', Adx(:)', 'gradf', gradf(:)', 'f', f, 'tau', ...
                          tau, 'epsilon', epsilon, 'alpha', alpha, ...
                          'beta', beta, 's', s, 'backiter', BI), jopts);
    end
    
    % set up for next iteration
    x = xp; u = up; r = rp; fu1 = fu1p; fu2 = fu2p; fe = fep; f = fp;
    
    lambda2 = -(gradf'*[dx; du]);
    stepsize = s*norm([dx; du]);
    
    if opts.verbose
      fprintf(['Newton iter = %d, Functional = %8.3f, '
                    'Newton decrement = %8.3f, Stepsize = %8.3e'], ...
                   niter, f, lambda2/2, stepsize);
      
      fprintf('                CG Res = %8.3e, CG Iter = %d', cgres, cgiter);
    end
    
    if (lambda2/2 < newtontol)
      return
    end
    
  end
end

function [xp, up, rp, fu1p, fu2p, fep, fp, backiter] = line_search(x, u, r, dx, du, Adx, gradf, f, tau, ...
                                    epsilon, alpha, beta, s)
% backtracking line search
  for backiter=1:32
    xp = x + s*dx;  
    up = u + s*du;
    rp = r + s*Adx;
    
    fu1p = xp - up;  
    fu2p = -xp - up;  
    fep = 1/2*(rp'*rp - epsilon^2);
    
    fp = sum(up) - (1/tau)*(sum(log(-fu1p)) + sum(log(-fu2p)) + log(-fep));
    flin = f + alpha*s*(gradf'*[dx; du]);
    
    s = beta*s;
    suffdec = (fp <= flin);
    if suffdec
      return
    end
  end
  fprintf(['Backtracking line search failed, returning previous ' ...
           'iterate.  (See Section 4 of notes for more information.)\n']);
  fprintf('fp = %f, flin = %f\n', fp, flin)
  keyboard
  xp = x;  
  up = u;
  rp = r;
  
end

function s = find_max_step(dx, du, Adx, fu1, fu2, r, epsilon)

% maximum step size that stays in the interior
  ifu1 = find((dx-du) > 0); 
  ifu2 = find((-dx-du) > 0);
  aqe = Adx'*Adx;   
  bqe = 2*r'*Adx;   
  cqe = r'*r - epsilon^2;
  smax = min(1,min([...
    -fu1(ifu1)./(dx(ifu1)-du(ifu1)); 
    -fu2(ifu2)./(-dx(ifu2)-du(ifu2)); ...
    (-bqe+sqrt(bqe^2-4*aqe*cqe))/(2*aqe)
                   ]));
  s = (0.99)*smax;
  
end



function [dx, du, gradf, cgres, gd] = compute_descent(fu1, fu2, r, fe, ...
                                                  tau,cgtol, cgmaxiter,A, At)
 
  % This is solving a sytem of (block)equations before passing to
  % cgsolve. Ie, the entire H_z of eq. (8) does not go in there:
  % [ H_xx, H_xu] [dx] = [delFx]
  % [ H_ux, H_uu] [du]   [delFu]
  % Then
  % du = Huu^-1 delFu - Huu^-1Hux*dx
  % and cgsolve solves the system of euqations for dx:
  % 
  % (Hxx - HxuHuu^-1Hux)x = delFu - HxuHuu^-1 *deltFu
  % The ananymous function is the LHS and w1p is the RHS.
  atr = At(r); % atr = A'*r
  ntgz = 1./fu1 - 1./fu2 + 1/fe*atr;
  ntgu = -tau - 1./fu1 - 1./fu2;
  gradf = -(1/tau)*[ntgz; ntgu];
  
  sig11 = 1./fu1.^2 + 1./fu2.^2;
  sig12 = -1./fu1.^2 + 1./fu2.^2;
  sigx = sig11 - sig12.^2./sig11;
    
  w1p = ntgz - sig12./sig11.*ntgu;
%   if (largescale)
    h11pfun = @(z) sigx.*z - (1/fe)*At(A(z)) + 1/fe^2*(atr'*z)*atr;
    [dx, cgres, cgiter] = cgsolve(h11pfun, w1p, cgtol, cgmaxiter, 0);
    du = (1./sig11).*ntgu - (sig12./sig11).*dx;    
  
   gd = struct('atr', atr, 'ntgu', ntgu, 'sig11', sig11, 'sig12', sig12,...
     'w1p', w1p, 'sigx', sigx); 
end
