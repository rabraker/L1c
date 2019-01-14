function build_dct_test_data(test_data_root)
  % Create a small set of test data
  jopts.FloatFormat = '%.20f';
  
  EMx_path = fullfile(test_data_root, 'dct_small_EMx.json');
  MtEty_path = fullfile(test_data_root, 'dct_small_MtEty.json');
  MtEt_Emx_path = fullfile(test_data_root, 'dct_small_MtEt_EMx.json');
  N = 50;
  Ts = 1/50;
  To = 0.25;
  
  omega =2*pi/To;
  
  k = [ 0:N-1]';
  
  x = sin(omega * Ts *k);
  
  pix_idx = [2, 10, 15, 20, 25, 30, 35, 40, 45, 49];
  
  yy = L1qcTestData.Afun_dct(x, pix_idx);
  
  jopts.FileName = EMx_path;
  data = struct('x0', x(:)', 'x1', yy(:)', 'pix_idx', pix_idx(:)'-1);
  savejson('', data, jopts);
  
  x_AtA = L1qcTestData.Atfun_dct(yy, pix_idx, N);

  
  
  jopts.FileName = MtEty_path;
  data = struct('x0', x_AtA(:)', 'x1', yy(:)', 'pix_idx', pix_idx(:)'-1);
  savejson('', data, jopts);
  
  jopts.FileName = MtEt_Emx_path;
  data = struct('x0', x(:)', 'x1', x_AtA(:)', 'pix_idx', pix_idx(:)'-1);
  savejson('', data, jopts);
  
  
  %%
  % Another small example with randomly generated data.

  Nx = 50;
      
  pix_mask_vec = ones(Nx,1);
  pix_idx = find(pix_mask_vec > 0.5);
  Ny = length(pix_idx);
  
  % y = E*M*x
  img_vec = randn(Nx, 1);
  
  % y, set of measurements. have to remove all the spots we didn't sample.
  y_vec = img_vec(pix_idx);
  A = @(x)idct(x);
  At = @(x)dct(x);
  
  EMx = A(img_vec);
  MtEty = At(y_vec);
  MtEt_EMx = At(A(img_vec));
  
  jopts.FloatFormat = '%.20f';
  jopts.FileName = fullfile(test_data_root, 'dct_small_rand.json');
  savejson('', struct('x_in', img_vec(:)', 'y_in', y_vec(:)',...
    'EMx', EMx(:)', 'MtEty', MtEty(:)', 'MtEt_EMx', MtEt_EMx(:)', 'pix_idx', pix_idx(:)'-1), jopts);
  
  
  %%  
  % Now, Use data from an actual compressed sensing situation. 
  
  img_dat = load('test_image_data.mat');
  xorig = img_dat.xorig;
  pix_idx = img_dat.pix_idx;
  N = length(xorig);
  
  A = @(x) L1qcTestData.Afun_dct(x, pix_idx); % E*M
  At = @(x) L1qcTestData.Atfun_dct(x, pix_idx, N); %E^T*M^T
  
  b = xorig(pix_idx);
  
  EMx = A(xorig);
  MtEty = At(b);
  
  MtEt_EMx = At(A(xorig));
  
  jopts.FloatFormat = '%.20f';
  jopts.FileName = fullfile(test_data_root, 'dct_large.json');
  savejson('', struct('x_in', xorig(:)', 'y_in', b(:)',...
    'EMx', EMx(:)', 'MtEty', MtEty(:)', 'MtEt_EMx', MtEt_EMx(:)', 'pix_idx', pix_idx(:)'-1), jopts);
  

end