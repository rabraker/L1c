function [npass, nfail, nskip] = test_mex_interface(skip_all, data_dir)
% Test cases for the mex bindings. data_dir should hold the test data and in
% particular 'example_img_data.json'
% 
  
  fpath = fullfile(data_dir, 'example_img_data.json');  
  
  cases = {...
    @()check_breg_TV_mex(fpath),...
    @()check_nesta_dctTV_notv(fpath),...
    @()check_nesta_dctTV(fpath),...
    @()check_l1qc_dct_mex(fpath),...
           };

  if ~skip_all
    [npass, nfail, nskip] = L1cMexTesting.run_suite(cases);
  else
    npass = 0;
    nfail = 0;
    nskip = length(cases);
  end

end

function check_nesta_dctTV(fpath)
  verbose = 0;
  fid = fopen(fpath, 'r');
  dat_json = fscanf(fid, '%s');
  fclose(fid);
  
  dat = jsondecode(dat_json);
  
  opts = nesta_opts('alpha_v', 0.5, 'alpha_h', .2, 'verbose', 0, 'dct_mode', 1, ...
                    'bp_mode', 'analysis');
  T = .20; %sparsity percent.
  if dat.one_based_index == 1
    pix_idx = dat.pix_idx;
  else
    pix_idx = dat.pix_idx+1;
  end
  % For the test below, we must sparsify the signal.
  x_orig = dat.x_orig;
  
  [x_est, status] = nesta_dctTV(dat.mrow, dat.mcol, x_orig(pix_idx), pix_idx, opts);

  L1cMexTesting.assert_eq(status, 0);

  % assert_bp_properties(x_og_sparse, x_est(:), pix_idx, T_idx, TC_idx, opts.sigma)
  % TODO: what kind of properties hold in the analysis case? This passes, but
  % I dont know why.
end


function check_nesta_dctTV_notv(fpath)
  verbose = 0;
  fid = fopen(fpath, 'r');
  dat_json = fscanf(fid, '%s');
  fclose(fid);
  
  dat = jsondecode(dat_json);
  
  opts = nesta_opts('alpha_v', 0, 'alpha_h', 0, 'verbose', 0, 'dct_mode', 1);
  T = .20; %sparsity percent.
  if dat.one_based_index == 1
    pix_idx = dat.pix_idx;
  else
    pix_idx = dat.pix_idx+1;
  end
  % For the test below, we must sparsify the signal.
  x_orig = dat.x_orig;
  [z_orig, T_idx, TC_idx] = sparsify(x_orig, T);
  x_og_sparse = idct(z_orig);
  
  [x_est, status] = nesta_dctTV(dat.mrow, dat.mcol, x_og_sparse(pix_idx), pix_idx, opts);

  L1cMexTesting.assert_eq(status, 0);

  assert_bp_properties(x_og_sparse, x_est(:), pix_idx, T_idx, TC_idx, opts.sigma)
  % TODO: implement the property check we use in the c-testsuite.
end


function check_l1qc_dct_mex(fpath)
  verbose = 0;
  fid = fopen(fpath, 'r');
  dat_json = fscanf(fid, '%s');
  fclose(fid);
  
  dat = jsondecode(dat_json);
  T = .20; %sparsity percent.
  epsilon = 0.1;
  opts = l1qc_dct_opts('verbose', verbose, 'l1_tol', 1e-5,...
                       'epsilon', epsilon, 'mu', 10);

  if dat.one_based_index == 1
    pix_idx = dat.pix_idx;
  else
    pix_idx = dat.pix_idx+1;
  end
  % For the test below, we must sparsify the signal.
  x_orig = dat.x_orig;
  [z_orig, T_idx, TC_idx] = sparsify(x_orig, T);
  x_og_sparse = idct(z_orig);
  
  [x_est, LBRes]= l1qc_dct(dat.mtot, 1, x_og_sparse(pix_idx), pix_idx, opts);
  L1cMexTesting.assert_eq(LBRes.status, 0);
  
  assert_bp_properties(x_og_sparse, x_est, pix_idx, T_idx, TC_idx, epsilon)
  
end

function assert_bp_properties(x_orig, x_est, pix_idx, T_idx, TC_idx, sigma)
% From Stable Signal Recovery from Incomplete and Inaccurate Measurements,
% Candes, Romberg, Tao 2005, it looks like we should be able to verify
% (cir. (9)-(10)):
%    a. ||x_opt||_1 <= ||x_act||_1
%    b. ||Ax_opt - Ax_act||_1 <= 2*\epsilon
%    c. Let h = x_opt - x_act. Let h_T = h[idx_supp], and
%    h_TC = h[(1:m)!=idx_supp]. Then
%    ||h_TC||_1 <= ||h_T||_1
%
%    d. From their numerical experiments, they suggest that
%       ||x_opt - x_act|| < C*eps, with C<2 (see also eq. (5)).
%       with eps^2 = sigma^2(n + lambda*sqrt(2*n)). Though table 1
%       shows this, I dont get the same result repeating that
%       experiment, even with their software. It seems that C~=3.5
%
  
  z_orig = dct(x_orig);
  z_est = dct(x_est);  
  
  % a. ||x_opt||_1 <= ||x_act||_1
    nrm_z_orig =  norm(z_orig, 1);
    nrm_z_opt = norm(z_est, 1);

    L1cMexTesting.assert_lt(nrm_z_opt, nrm_z_orig);

  % b). ||Az_opt - Az_act||_1 <= 2*\sigma
  %   = ||E*x_opt - E*x_act||

  y_err = x_orig(pix_idx) - x_est(pix_idx);
  L1cMexTesting.assert_lt(norm(y_err), 2*sigma);
  
  % c) c. Let h = x_opt - x_act. Let h_T = h[idx_supp], and
  %     h_TC = h[(1:m)!=idx_supp]. Then
  %     ||h_TC||_1 <= ||h_T||_1
  h = z_orig - z_est;
  h_T = h(T_idx);
  h_TC = h(TC_idx);
  L1cMexTesting.assert_lt(norm(h_TC, 1), norm(h_T, 1));
    
end

function check_breg_TV_mex(fpath)
  fid = fopen(fpath, 'r');
  dat_json = fscanf(fid, '%s');
  fclose(fid);
  
  dat = jsondecode(dat_json);
  
  n = sqrt(dat.mtot);
  m = n;
  x = dat.x_orig(:) + rand(dat.mtot, 1);
  img = reshape(x, n, n);
  mu = 5;
  tol = 0.001;
  
  Xopt = breg_anistropic_TV(img, mu, tol, 1000);
  
  L1cMexTesting.assert_eq(size(Xopt, 1), n);
  L1cMexTesting.assert_eq(size(Xopt, 2), m);
  
  err = norm(Xopt(:) - x(:))/norm(x(:));
  L1cMexTesting.assert_lt(err, 1.0);
end

function [Z, T_idx, TC_idx] = sparsify(x, T)
  [N, M] = size(x);
  img_vec = reshape(x', N*M,1);
  k = floor(T*N*M);
  Z = dct(img_vec);
  [z_srt] = sort(abs(Z), 'descend');
  z_thresh = z_srt(k);
  Z(abs(Z) < z_thresh) = 0;
  
  T_idx = find(Z~=0);
  TC_idx = find(Z==0);
  
end
