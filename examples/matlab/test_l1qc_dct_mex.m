% Helper function for autotools testsuite.
function [tm_mex, x_est, LBRes] = test_l1qc_dct_mex(fpath, verbose)
  fid = fopen(fpath, 'r');
  dat_json = fscanf(fid, '%s');
  fclose(fid);
  
  dat = jsondecode(dat_json);
  
  opts = l1qc_dct_opts('verbose', verbose, 'l1_tol', 1e-5,...
                       'epsilon', 0.1, 'mu', 10);
  tic
  if dat.one_based_index == 1
    pix_idx = dat.pix_idx;
  else
    pix_idx = dat.pix_idx+1;
  end
  [x_est, LBRes]= l1qc_dct(dat.N, 1, dat.b, pix_idx, opts);
  tm_mex = toc;
  fprintf('mex file time: %.4f\n', tm_mex);

  if LBRes.status > 0
    exit(LBRes.status)
  end

end
