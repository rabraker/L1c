function build_cgsolve_test_data(test_data_root)
% Build and save test data for cgsolve.c
% 
% The JSON file will be saved to test_data_root/''cgsolve_small01.json'
%
  
  
  cg_test_data_small_path = fullfile(test_data_root, 'cgsolve_small01.json');
  tol = 1e-6;
  maxiter = 100;
  verbose = 1;
  
  % A slightly larger example.
  rng(2);
  N = 50;
  A_mat = rand(N,N);
  A_mat = A_mat + A_mat' + 6*eye(N);
  
  assert(min(eig(A_mat)) > 0)
  
  b = rand(N,1);
  A_fun = @(x) A_mat *x;
  x = L1qcTestData.cgsolve(A_fun, b, tol, maxiter, verbose);
  
  
  A_row = [];
  for i=1:size(A_mat, 1)
    A_row = [A_row, A_mat(i,:)];
  end
  
  data = struct('A', A_row, 'b', b', 'x', x', 'tol', tol, 'max_iter', maxiter);
  
  savejson('', data, cg_test_data_small_path);

end

