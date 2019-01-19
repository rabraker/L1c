% Builds test data for the prelude of l1qc_logbarrier.m
%


function xp = build_logbarrier_test_data(data_root, lbiter)  
  img_dat = load(fullfile(data_root, 'test_image_data.mat'));
  xorig = img_dat.xorig;
  pix_idx = img_dat.pix_idx;
  
  A = @(x) L1qcTestData.Afun_dct(x,pix_idx); % E*M
  At = @(x) L1qcTestData.Atfun_dct(x,pix_idx, length(xorig)); %E^T*M^T

  
  b = xorig(pix_idx);
  x0 = At(b);
 
  epsilon = 0.1;
  lbtol = 1e-3;
  mu = 10;
  cgtol = 1e-9;
  cgmaxiter = 200;

  newtontol = lbtol;
  newtonmaxiter = 50;
  verbose = 1;
  N = length(x0);

  % starting point --- make sure that it is feasible
  if (norm(A(x0)-b) > epsilon)
    error('Starting point infeasible; using x0 = At*inv(AAt)*y.');
  end
  
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
  x_k = x0;
  for ii = 1:lbiter
   
    [xp, up, ntiter] = L1qcTestData.l1qc_newton(x_k, u, A, At, b, epsilon, tau,...
      newtontol, newtonmaxiter, cgtol, cgmaxiter, ii, verbose, 0);
    totaliter = totaliter + ntiter;
    
    fprintf(['\nLog barrier iter = %d, l1 = %.3f, functional = %8.3f, ',...
      'tau = %8.3e, total newton iter = %d\n'], ...
            ii, sum(abs(xp)), sum(up), tau, totaliter);
    
    fname_iter = sprintf('lb_test_data_iter_%d.json', ii);
    jopts.FileName = fullfile(data_root, fname_iter);
    savejson('', struct('x0', x0(:)', 'b', b(:)',...
      'epsilon', epsilon, 'tau', tau, 'lbiter', lbiter, 'lbtol', lbtol,...
      'newtontol', newtontol, 'newtonmaxiter', newtonmaxiter,...
      'cgtol', cgtol, 'cgmaxiter', cgmaxiter, 'mu', mu, 'pix_idx', pix_idx(:)'-1,...
      'xp', xp(:)', 'up', up(:)', 'tau_next', mu*tau), jopts);
    
    tau = mu*tau;
    x_k = xp;
    u = up;
    
  end

end
