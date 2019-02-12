% This will generate an opts struct suitible for passing to l1qc(). 
function [tm_mex, x_est, LBRes] = test_l1qc_dct_mex(fpath, verbose)
  
  dat=loadjson(fpath);
  opts = l1qc_dct_opts('verbose', verbose);
  tic
  [x_est, LBRes]= l1qc_dct_mex(dat.N, dat.b, dat.pix_idx-1, opts);
  tm_mex = toc;
  fprintf('mex file time: %.4f\n', tm_mex);
end
