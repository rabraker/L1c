function build_newton_init_data(data_root)
  img_dat = load(fullfile(test_data_root, 'test_image_data.mat'));


  xorig = img_dat.xorig;
  pix_idx = img_dat.pix_idx;
  
  A = @(x) L1qcTestData.Afun_dct(x,pix_idx); % E*M
  At = @(x) L1qcTestData.Atfun_dct(x,pix_idx, length(xorig)); %E^T*M^T

  
  b = xorig(pix_idx);
  x0 = At(b);
  
  
  epsilon = 0.1;
      
  lbtol = 1e-3;
  mu = 10;
  N = length(x0);
  
      % starting point --- make sure that it is feasible
  if (norm(A(x0)-b) > epsilon)
    error('Starting point infeasible; using x0 = At*inv(AAt)*y.');
  end
  
  [u, lbiter, tau] = build_init_dat(x0, mu, lbtol);

  fprintf('Original l1 norm = %.3f, original functional = %.3f\n', sum(abs(x0)), sum(u));


  jopts.FileName = fullfile(data_root, 'newton_init_data.json');
  jopts.FloatFormat = '%.20f';
  
  savejson('', struct('x', x0(:)', 'u', u(:)', 'b', b(:)',...
    'epsilon', epsilon, 'tau', tau, 'lbiter', lbiter, 'lbtol',...
    lbtol, 'mu', mu, 'pix_idx', pix_idx(:)'-1), jopts);
  
  % ------------------------------------------------------------------
  % Add regression data to test for the bug where I had
  % max(xmax, x0[0]), rather than x0[i]
end



function [u, lbiter, tau] = build_init_dat(x0, mu, lbtol)
  
  N = length(x0);
  u = (0.95)*abs(x0) + (0.10)*max(abs(x0));
  
  % choose initial value of tau so that the duality gap after the first
  % step will be about the origial norm
  tau = max((2*N+1)/sum(abs(x0)), 1);

  lbiter = ceil((log(2*N+1)-log(lbtol)-log(tau))/log(mu));
  fprintf('Number of log barrier iterations = %d\n\n', lbiter);
  
  
end