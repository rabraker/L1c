% This will generate an opts struct suitible for passing to l1qc(). 
function [tm_mex, x_est, LBRes] = test_l1qc_dct_mex(fpath, verbose)
  
  dat=loadjson(fpath);
  opts = l1qc_dct_opts('verbose', verbose);
  tic
  if dat.one_based_index == 1
    pix_idx = dat.pix_idx;
  else
    pix_idx = dat.pix_idx+1;
  end
  [x_est, LBRes]= l1qc_dct_mex(dat.N, 1, dat.b, pix_idx, opts);
  tm_mex = toc;
  fprintf('mex file time: %.4f\n', tm_mex);

  if LBRes.status > 0
    exit(LBRes.status)
  end

end
