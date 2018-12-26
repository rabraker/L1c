function build_newton_init_data(data_root)
  
  fpath = '/home/arnold/matlab/afm-cs/matlab-code/notes/data/cs_sim_CS20NG.mat';
  addpath ~/matlab/afm-cs/matlab-code/functions
  addpath ~/matlab/afm-cs/reconstruction/BP
  
  dat = load(fpath);
  cs_sim = dat.cs_sim;
  pix_mask_vec = PixelMatrixToVector(cs_sim.pix_mask);
  pix_idx = find(pix_mask_vec>0.5);
  
  b = PixelMatrixToVector(cs_sim.Img_original);
  b = b(pix_idx); 
  b = b/max(b);
  
  A = @(x) L1qcTestData.IDCTfun(x,pix_mask_vec);
  At = @(x) L1qcTestData.DCTfun(x,pix_mask_vec);
  x0 = round(At(b), 15);
  
  epsilon = 0.1;
      
  lbtol = 1e-3;
  mu = 10;
  N = length(x0);
  
      % starting point --- make sure that it is feasible
  if (norm(A(x0)-b) > epsilon)
    error('Starting point infeasible; using x0 = At*inv(AAt)*y.');
  end
  
  [u, lbiter, tau] = build_init_dat(x0)

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



function [u, lbiter, tau] = build_init_dat(x0)
  
  u = (0.95)*abs(x0) + (0.10)*max(abs(x0));
  
  % choose initial value of tau so that the duality gap after the first
  % step will be about the origial norm
  tau = max((2*N+1)/sum(abs(x0)), 1);

  lbiter = ceil((log(2*N+1)-log(lbtol)-log(tau))/log(mu));
  fprintf('Number of log barrier iterations = %d\n\n', lbiter);
  
  
end