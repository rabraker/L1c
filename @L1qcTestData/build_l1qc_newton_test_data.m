function build_l1qc_newton_test_data(test_data_root)
%   fpath = '/home/arnold/matlab/afm-cs/matlab-code/notes/data/cs_sim_CS20NG.mat';
%   addpath ~/matlab/afm-cs/matlab-code/functions
% 
%   dat = load(fpath);
%   cs_sim = dat.cs_sim;
% 
% 
%   pix_mask_vec = L1qcTestData.pixmat2vec(cs_sim.pix_mask);
% 
%   % y = E*M*x
%   y_vec = L1qcTestData.pixmat2vec(cs_sim.Img_sub_sampled);
%   % y, set of measurements. have to remove all the spots we didn't sample.
%   y_vec = y_vec(pix_mask_vec>0.5);
%   y_vec = y_vec/max(y_vec); % normalize to 1
%   pix_idx = find(pix_mask_vec > 0.5); % for saving, not computation. 
  img_dat = load(fullfile(test_data_root, 'test_image_data.mat'));
  xorig = img_dat.xorig;
  pix_idx = img_dat.pix_idx;
  
  A = @(x) L1qcTestData.Afun_dct(x,pix_idx); % E*M
  At = @(x) L1qcTestData.Atfun_dct(x,pix_idx, length(xorig)); %E^T*M^T

  
  b = xorig(pix_idx);
  x0 = At(b);
  

  opts.epsilon = 0.1;
  u = (0.95)*abs(x0) + (0.10)*max(abs(x0));
  N = length(x0);
  % choose initial value of tau so that the duality gap after the first
  % step will be about the origial norm
  opts.tau = max((2*N+1)/sum(abs(x0)), 1);
  opts.cgtol = 1e-8;
  opts.cgmaxiter = 200;
  opts.lbtol = 1e-3;
  opts.newtontol = opts.lbtol;
  opts.newtonmaxiter = 50;
  opts.verbose = 1;

  
  

  [xp2, up2, ntiter2] = l1qc_newton_local(x0, u, A, At, b, opts, test_data_root, pix_idx);
  
end


function [x, u, niter] = l1qc_newton_local(x0, u0, A, At, b, opts, data_root, pix_idx) 
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
  [r, fu1, fu2, fe, f] =  f_eval(A, x, u, b, epsilon, tau);
  
  jopts.FloatFormat ='%.20e';
    % Initial values for cgsolve init.
    
  if opts.verbose
    fprintf(['Newton-iter | Functional | Newton decrement |  Stepsize  ' ...
             '|  cg-res | cg-iter | backiter |  s0    |    s     |\n']);
  end    
  dx = 0*x0;
  cgit_tot = 0;
  s = 1;  

  for niter=1:newtonmaxiter
    if niter==3
      save_on = true;
    else
      save_on = false;
    end
    jopts.FileName = fullfile(data_root);
    descent_files = {fullfile(data_root, 'descent_data.json'),...
                     fullfile(data_root, 'hp11_fun_data.json')};
                   
    [dx, du, gradf, cgres, cgiter] = compute_descent(fu1, fu2, r, fe, tau,...
                                 cgtol, cgmaxiter, A, At, save_on, ...
                                                     jopts, descent_files, pix_idx, 0*dx);
   
    if (cgres > 1/2)
      disp('Cannot solve system.  Returning previous iterate.  (See Section 4 of notes for more information.)');
      xp = x;  up = u;
      return
    end
    
    Adx = A(dx);
    
    % ------------- get the maximum step size.------------------
    jopts.FileName = fullfile(data_root,'find_max_step_data.json');
    s = find_max_step(dx, du, Adx, fu1, fu2, r, epsilon, save_on, jopts);
    s0 = s;
    % ------------- Line search -------------------------------
     jopts.FileName = fullfile(data_root, 'line_search_data.json');
    % [xp, up, rp, fu1p, fu2p, fep, fp, BI]
    [x, u, r, fu1, fu2, fe, f, BI, s] = line_search(x, u, r, b, dx, du, A, gradf, f, tau, ...
                                    pix_idx, epsilon, alpha, beta, s, save_on, ...
                                                      jopts);
    
    % set up for next iteration
    % x = xp; u = up; r = rp; fu1 = fu1p; fu2 = fu2p; fe = fep; f = fp;
    
    lambda2 = -(gradf'*[dx; du]);
    stepsize = s*norm([dx; du]);
    
    if opts.verbose
      % fprintf('Newton iter |  Functional | Newton decrement | Stepsize | cg-res | backiter|  s \n');
      %            NI         fcnl         dec            sz     cgr       cgI        BI       s
      fprintf('     %3d       %8.3g       %08.3g       % 8.3e   %08.3g     %3d       %2d   %.3g    %.3g \n',...
        niter, f, lambda2/2, stepsize, cgres, cgiter, BI, s0,  s);
    end
    
    
    if (lambda2/2 < newtontol)
      return
    end
    
  end
end

function [xp, up, rp, fu1p, fu2p, fep, fp, backiter, s] = line_search(x, u, r, b, dx, du, A, gradf, f, tau, ...
                                    pix_idx, epsilon, alpha, beta, s_init, save_on, ...
                                                    jopts)
% backtracking line search
  n = length(x);
  s = s_init;
  for backiter=1:32
    xp = x + s*dx;  
    up = u + s*du;
    %rp = r + s*Adx;
%     rp = A(xp) - b;
%     fu1p = xp - up;  
%     fu2p = -xp - up;  
%     fep = 1/2*(rp'*rp - epsilon^2);
%     fp = sum(up) - (1/tau)*(sum(log(-fu1p)) + sum(log(-fu2p)) + log(-fep));
  
    [rp, fu1p, fu2p, fep, fp] =  f_eval(A, xp, up, b, epsilon, tau);
    flin = f + alpha*s*(gradf'*[dx; du]);
    
    flx = gradf(1:n)'*dx;
    flu = gradf(n+1:end)'*du;
    
    
    suffdec = (fp <= flin);
    if suffdec
      break
    end
    s = beta*s;
  end
  
  
  if save_on

    savejson('', struct('xp', xp(:)', 'up', up(:)', 'rp', rp(:)', 'fu1p', fu1p(:)', ...
                        'fu2p', fu2p(:)', 'fep', fep(:)', 'fp', fp(:)', 'x', ...
                        x(:)', 'u', u(:)', 'r', r(:)', 'b', b(:)', 'dx', dx(:)', 'du', du(:)', ...
                        'gradf', gradf(:)', 'f', f, 'tau', tau,...
                         'pix_idx', pix_idx(:)'-1,'epsilon', epsilon, 'alpha', alpha, ...
                        'beta', beta, 's', s, 's_init', s_init, 'backiter', backiter,...
                        'flx', flx, 'flu', flu, 'flin', flin), jopts);
  end
if suffdec
  return
else
    fprintf(['Backtracking line search failed, returning previous ' ...
           'iterate.  (See Section 4 of notes for more information.)\n']);
  fprintf('fp = %f, flin = %f\n', fp, flin)
  keyboard
  xp = x;  
  up = u;
  rp = r;  
end
  
end

function s = find_max_step(dx, du, Adx, fu1, fu2, r, epsilon, save_on, jopts)

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
  
  if save_on
    savejson('', struct('dx', dx(:)', 'du', du(:)', 'Adx', Adx(:)', 'fu1', fu1(:)', ...
                        'fu2', fu2(:)', 'r', r(:)', 'epsilon', epsilon, ...
                        'smax', s), jopts);
  end
  
end



function [dx, du, gradf, cgres, cgiter] = compute_descent(fu1, fu2, r, fe, ...
                         tau,cgtol, cgmaxiter,A, At, save_on, jopts, ...
                                                    path_spec, pix_idx, dx0)
 
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
  
  verbose = 0;
  if save_on
    verbose = 1;
  end
  h11pfun = @(z) sigx.*z +(-(1/fe)*At(A(z)) + 1/fe^2*(atr'*z)*atr);
  [dx, cgres, cgiter] = L1qcTestData.cgsolve(h11pfun, w1p, cgtol, ...
                                             cgmaxiter, verbose, dx0);
  du = (1./sig11).*ntgu - (sig12./sig11).*dx;

%    gd = struct('atr', atr, 'ntgu', ntgu, 'sig11', sig11, 'sig12', sig12,...
%      'w1p', w1p, 'sigx', sigx); 
%    
   if save_on
     jopts.FileName = path_spec{1};
     savejson('', struct('dx0', dx0(:)', 'dx', dx(:)', 'du', du(:)', 'gradf', gradf(:)', 'fu1', ...
                        fu1(:)', 'fu2', fu2(:)', 'r', r(:)', 'atr', atr(:)',...
                        'ntgu', ntgu(:)', 'sig11', sig11(:)', 'sig12', sig12(:)',...
                        'w1p', w1p(:)', 'sigx', sigx(:)', 'pix_idx', pix_idx(:)'-1,...
                        'fe', fe, 'tau', tau, 'cgtol', cgtol, 'cgmaxiter', cgmaxiter,...
                        'cgres', cgres, 'cgiter', cgiter), jopts);
                      
     % hp11_fun data
     jopts.FileName = path_spec{2};
     z_typ = w1p; % From line 36 & 28 of cgsolve
%      h11p = @(z) At(A(z));% sigx.*z - (1/fe)*+ 1/fe^2*(atr'*z)*atr;
     y_exp = h11pfun(z_typ);
     savejson('', struct('z', z_typ(:)', 'atr', atr(:)', 'fe', fe, 'sigx',...
       sigx(:)', 'y_exp', y_exp(:)', 'pix_idx', pix_idx(:)'-1), jopts);
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

