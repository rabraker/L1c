function build_cgsolve_test_data(test_data_root)
% Build and save test data for cgsolve.c
% 
% The JSON file will be saved to test_data_root/''cgsolve_small01.json'
%
  
  
  cg_test_data_small_path = fullfile(test_data_root, 'cgsolve_small01.json');
  tol = 1e-9;
  maxiter = 50;
  verbose = 1;
  
  % A slightly larger example.
  rng(2);
  N = 50;
  A_mat = rand(N,N)*50;
  A_mat = round(A_mat, 0);
  A_mat = A_mat + A_mat' + 600*eye(N);
  
  
  % NOTE: errors between matlab cgsolve and c cgsolve for this small problem
  % seem to come from errors between the two matrix-vector multiplications. I
  % tested this by stopping each solver in the first iteration immediatly after
  % the q = A*d step and printing out d and q to 16 decimal places (in the c
  % code). Copying and pasting that into matlab, I get that d_c - d_matlab = 0
  % everywhere, but q_c - q_matlab has an error of about 1e-13. I can only
  % guess that this is because the two dgemv function must do things in
  % slightly different order. I don't really know, but that's not what I'm
  % trying to test here. So, Im going to multiply the random data up and round
  % to zero decimal places. The resulting vector solution are now accurate to
  % within 
  assert(min(eig(A_mat)) > 0)
  
  b = rand(N,1)*50;
  b = round(b, 0);
  
  
  A_fun = @(x) A_mat *x;
  x = L1qcTestData.cgsolve(A_fun, b, tol, maxiter, verbose);
  
  
  A_row = [];
  for i=1:size(A_mat, 1)
    A_row = [A_row, A_mat(i,i:end)];
  end
  
  data = struct('A', A_row, 'b', b', 'x', x', 'tol', tol, 'max_iter', maxiter);
  
  jopts.FileName = cg_test_data_small_path;
  jopts.FloatFormat = '%.20f';
  savejson('', data, jopts);

  
  % ------------------------------------------------ %
  rng(1);
  N = 50;
  A = rand(N,N);
  
  A = (A + A')/2 + 6*eye(N);
  
  assert(min(eig(A)) >0)
  
  
  x = rand(N,1);
  
  y = A*x;
  
  Arow = [];
  for i=1:size(A, 1)
    Arow = [Arow, A(i,i:end)];
  end
  
  jopts.FileName = 'test_data/ax_sym.json';
  jopts.FloatFormat = '%.20f';
  
  savejson('', struct('A', Arow(:)', 'x', x(:)', 'y', y(:)'), jopts);
end

