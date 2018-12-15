function build_dct_test_data(test_data_root)
  % Create a small set of test data

  EMx_path = fullfile(test_data_root, 'dct_small_EMx.json');
  MtEt_path = fullfile(test_data_root, 'dct_small_MtEt_EMx.json');
  MtEt_Emx_path = fullfile(test_data_root, 'dct_small_MtEt_EMx.json');
  N = 50;
  Ts = 1/50;
  To = 0.25;
  
  omega =2*pi/To;
  
  k = [ 0:N-1]';
  
  x = sin(omega * Ts *k);
  
  pix_idx = [2, 10, 15, 20, 25, 30, 35, 40, 45, 49];
  
  Ny = length(pix_idx);
  pix_mask = zeros(N,1);
  pix_mask(pix_idx) = 1;
  
  yy = IDCTfun(x, pix_mask);
  
  data = struct('x0', x(:)', 'x1', yy(:)', 'pix_idx', pix_idx(:)'-1);
  savejson('', data, EMx_path);
  
  
  xx = DCTfun(yy, pix_mask);
  xx_save = xx;
   
  data = struct('x0', xx_save(:)', 'x1', yy(:)', 'pix_idx', pix_idx(:)'-1);
  savejson('', data, MtEt_path);
  %
  data = struct('x0', x(:)', 'x1', xx_save(:)', 'pix_idx', pix_idx(:)'-1);
  savejson('', data, MtEt_Emx_path);
  
end