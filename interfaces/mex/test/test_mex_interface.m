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
  
  opts = l1qc_dct_opts('verbose', verbose, 'l1_tol', 1e-5,...
                       'epsilon', 0.1, 'mu', 10);

  if dat.one_based_index == 1
    pix_idx = dat.pix_idx;
  else
    pix_idx = dat.pix_idx+1;
  end
  [x_est, LBRes]= l1qc_dct(dat.mtot, 1, dat.b, pix_idx, opts);

  L1cMexTesting.assert_eq(LBRes.status, 0);
  
  % TODO: implement the property check we use in the c-testsuite.
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

