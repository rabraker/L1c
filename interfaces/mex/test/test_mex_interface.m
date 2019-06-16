function pass_fail = test_mex_interface(data_dir)
% Test cases for the mex bindings. data_dir should hold the test data and in
% particular 'example_img_data.json'
% 
  fpath = fullfile(data_dir, 'example_img_data.json');

  MEX_RUN_CASE = getenv('MEX_RUN_CASE');

  cases = {...
    @()check_breg_TV_mex(fpath),...
    @()test_nesta_dctTV(fpath),...
    @()test_l1qc_dct_mex(fpath),...
           };

  pass_fail = L1cMexTesting.run_suite(cases);

end


function test_nesta_dctTV(fpath)
  verbose = 0;
  fid = fopen(fpath, 'r');
  dat_json = fscanf(fid, '%s');
  fclose(fid);
  
  dat = jsondecode(dat_json);
  
  opts = nesta_opts('alpha_v', 1, 'alpha_h', 1, 'verbose', 0);

  if dat.one_based_index == 1
    pix_idx = dat.pix_idx;
  else
    pix_idx = dat.pix_idx+1;
  end

  [x_est, status] = nesta_dctTV(dat.mrow, dat.mcol, dat.b, pix_idx, opts);

  L1cMexTesting.assert_eq(status, 0);

  % TODO: implement the property check we use in the c-testsuite.
end


function test_l1qc_dct_mex(fpath)
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
  
  % a. ||x_opt||_1 <= ||x_act||_1
  z_est = dct(x_est);
  nrm_z_opt = norm(z_est, 1);
  nrm_z_orig =  norm(z_orig, 1);
  L1cMexTesting.assert_lt(nrm_z_opt, nrm_z_orig);

  % b). ||Az_opt - Az_act||_1 <= 2*\epsilon
  %   = ||E*x_opt - E*x_act||
  y_err = x_og_sparse(pix_idx) - x_est(pix_idx);
  L1cMexTesting.assert_lt(norm(y_err), 2*epsilon);
  
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
