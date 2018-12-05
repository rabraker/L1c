A = [4, 1; 1 3];
b = [1;2];

tol = 1e-6;
maxiter = 100;

verbose = 1;

cgsolve_tmp(A, b, tol, maxiter, verbose)

%%

rng(2);
N = 50;
A = rand(N,N);
A = A + A' + 6*eye(N);

min(eig(A))

b = rand(N,1)
x=cgsolve_tmp(A, b, tol, maxiter, verbose)


A_row = [];
for i=1:size(A, 1);
  A_row = [A_row, A(i,:)];
end

data = struct('A', A_row, 'b', b', 'x', x', 'tol', tol, 'max_iter', maxiter);

savejson('', data, 'test_data/cgsolve_small01.json');